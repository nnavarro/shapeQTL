####################################### 
# shapeQTL mapping experiment with R
#
# Nicolas Navarro - 2013-2014
########################################  
multiv.scanone <- function(cross,pheno,mod.red,covar,back.qtl=NULL,test="Pillai"){
  # mod.red may be either formula or an externally fitted null model
  chr <- names(cross$geno)
  if ("X"%in%chr) warning("in case of the X chromosome, mapping may be wrong")
  #eventually re-fit the null model y ~ mu + covar
  if (class(mod.red)[1]=="formula"){
    fm.red <- as.formula(paste("pheno",deparse(mod.red[-2],width.cutoff=500L)))
    if (length(fm.red)!=3) stop("model too long: update covariate names")	
    mod.red <- lm(fm.red)
  }
  else fm.red <- formula(mod.red)
  SSCPerr.red <- crossprod(mod.red$residuals)
  rank.S <- qr(SSCPerr.red)$rank #min(n,2k-4) 
  if(is.null(back.qtl)) fm.full <- as.formula(paste(paste("pheno",deparse(fm.red[-2])),"qtl",sep="+")) 
  else fm.full <- as.formula(paste(paste("pheno",deparse(fm.red[-2])),"back.qtl", "qtl",sep="+"))
  
  result <- NULL
  if (class(cross)[1]=="bc"){
    for (j in chr){
      pr <- cross$geno[[j]]$prob
      map <- attr(pr,"map")
      pr <- pr[,,-dim(pr)[3],drop=TRUE]
      lod <- apply(pr,2,lm.shape.test, pheno, covar, fm.full, SSCPerr.red, mod.red$rank, rank.S, back.qtl, test)
      z <- data.frame(chr=rep(j,length(map)),pos=map,lod=lod)
      rownames(z) <- names(map)
      class(z) <- c("scanone","data.frame")
      result <- rbind(result,z)
    } 
    class(result) <- c("scanone","data.frame")
  }
  else if(class(cross)[1]=="f2"){
    #for intercross there are 3 test ({full or add} vs null and full vs add)
    for (j in chr){
      pr <- cross$geno[[j]]$prob
      map <- attr(pr,"map")
      if(is.null(back.qtl)) fm.add <- as.formula(paste(paste("pheno",deparse(fm.red[-2])),"Exp.A",sep="+"))
      else fm.add <- as.formula(paste(paste("pheno",deparse(fm.red[-2])),"back.qtl","Exp.A",sep="+"))
      lod.dom <- apply(pr,2,lm.shape.test.partial, pheno, covar, fm.full, fm.add, class(cross)[1],back.qtl)
      #if f2: additive model: Expected 0, 1 or 2 alleles of type B
      Exp.A <- pr[,,2]+2*pr[,,3]
      lod.add <- apply(Exp.A,2,lm.shape.test, pheno, covar, fm.full, SSCPerr.red, mod.red$rank, rank.S, back.qtl, test)
      pr <- pr[,,-dim(pr)[3],drop=TRUE]
      lod.full <- apply(pr,2,lm.shape.test, pheno, covar, fm.full, SSCPerr.red, mod.red$rank, rank.S, back.qtl, test)
      z <- data.frame(chr=rep(j,length(map)),pos=map, lod=lod.full,lod.add=lod.add,dom=lod.dom)
      rownames(z) <- names(map)
      class(z) <- c("scanone","data.frame")
      result <- rbind(result,z)
    }	
  }
  else {print("this kind of cross is not yet implemented")}
  return(result)
}
#-------
lm.shape.test <- function(qtl,pheno,covar,fm.full,SSCPerr.red,mod.red.rank,rank.E,back.qtl=NULL,test="Pillai"){
  #-> Quicker to call directly the C warper of the Fortran code and simplier than calling Fortran code directly 
  #1st Set design matrix x = [1,sexM,momF1,log.CS,a]
  x <- model.matrix(as.formula(deparse(fm.full[-2])))
  #At R version 3.1.1 a supplementary argument appears in the call...
  mod.full <- .Call(stats:::C_Cdqrls,x=x,y=pheno,tol=1e-07,FALSE)
  n.ind <- nrow(pheno)
  #qtl model
  dfeff <- mod.full$rank - mod.red.rank
  dferr <- n.ind - mod.full$rank
  SSCPerr.full <- crossprod(mod.full$residuals)	
  #partial F-test: Full model vs Reduced model
  SSCPfull <- SSCPerr.red - SSCPerr.full
  if (pmatch(test,"Pillai",nomatch=0)) return(Pillai.test(SSCPfull,SSCPerr.full,dfeff,dferr,rank.E))
  if (pmatch(test,"Lik.ratio",nomatch=0)) return(LikelihoodRatio.test(SSCPerr.full,SSCPerr.red,dfeff,dferr,rank.E))
  if (pmatch(test,"Hotelling.Lawley",nomatch=0)) return(Hotelling.test(SSCPfull,SSCPerr.full,dfeff,dferr,rank.E))	
}
#-------
lm.shape.test.partial <- function(qtl,pheno,covar,fm.full,fm.add,cross.type,back.qtl=NULL,test="Pillai"){
  #if f2: additive model: Expected 0, 1 or 2 alleles of type B
  if(cross.type =="f2") {
    if (is.na(dim(qtl)[3])) Exp.A <- qtl[,2]+2*qtl[,3]
    else Exp.A <- qtl[,,2]+2*qtl[,,3]
  }
  
  x <- model.matrix(as.formula(deparse(fm.full[-2])))
  mod.full <- .Call(stats:::C_Cdqrls,x=x,y=pheno,tol=1e-07,FALSE)
  
  x <- model.matrix(as.formula(deparse(fm.add[-2])))
  mod.add <- .Call(stats:::C_Cdqrls,x=x,y=pheno,tol=1e-07,FALSE)
  
  SSCPerr.red <- crossprod(mod.add$residuals)
  mod.red.rank <- mod.add$rank
  rank.E <- qr(SSCPerr.red)$rank
  #qtl model
  n.ind <- nrow(pheno)
  dfeff <- mod.full$rank - mod.red.rank
  dferr <- n.ind - mod.full$rank	
  SSCPerr.full <- crossprod(mod.full$residuals)	
  #partial F-test: Full model vs Reduced model
  SSCPfull <- SSCPerr.red - SSCPerr.full
  
  if (pmatch(test,"Pillai",nomatch=0)) return(Pillai.test(SSCPfull,SSCPerr.full,dfeff,dferr,rank.E))
  if (pmatch(test,"Lik.ratio",nomatch=0)) return(LikelihoodRatio.test(SSCPerr.full,SSCPerr.red,dfeff,dferr,rank.E))
  if (pmatch(test,"Hotelling.Lawley",nomatch=0)) return(Hotelling.test(SSCPfull,SSCPerr.full,dfeff,dferr,rank.E))	
}