#' fitqtl.shape is a function that allows estimating qtl effect.
#' @title fitqtl.shape
#' @author Nicolas Navarro
#' @param cross A cross object.
#' @return The function returns the qtl shape effects, eventually projected back in the tangent space.
#' @references xxxx
#' @seealso  \code{\link{fitqtl}}
#' @keywords shape effect
#' @note If cross is ....
#' @examples
#' eff <- fitqtl.shape(cross, pheno.col=1:2, qtl, covar=NULL, "y~covar")
#' @export

fitqtl.shape <- function(cross, pheno.col, qtl, covar=NULL, formula, 
                         method="hk", model="mvnorm",dropone=FALSE, get.ests=TRUE, run.checks=TRUE, forceXcovar=FALSE, pca=NULL) 
{
  if(!any(class(cross)=="cross")) stop("Cross should have class \"cross\".")
  if(!any(class(cross)%in%c("bc","f2"))) stop("Cross should be an \"f2\" or a \"bc\".")
  if(!any(class(cross)=="shape")) stop("Cross should have class \"shape\".")
  
  #qtl may be q qtl object or a summary from scanone.shape
  if(!any(class(qtl)%in%c("qtl","summary.scanone"))) stop("qtl argument should have class \"qtl\".")
  if(any(class(qtl)%in%"summary.scanone")){
    qtls.name <- find.pseudomarker(cross,chr=qtls$chr,pos=qtls$pos,where="prob")
    Q <- makeqtl(cross,qtls$chr,qtls$pos,qtls.name,what="prob")
  }else Q <- qtl
  n.gen <- Q$n.gen[1]
  if (any(Q$n.gen!=n.gen)) stop("QTLs have different number of alleles")
  
  if(is.character(pheno.col)) {
    if(length(which(colnames(cross$pheno)%in%pheno.col)))
      stop("Phenotype \"", pheno.col, "\" don't match")
    pheno.col <- which(colnames(cross$pheno)%in%pheno.col)
  }
  if (missing(formula)) fm.full <- paste(fm.red,paste("Q",1:Q$n.qtl,sep="",collapse="+"),sep="+")
  else fm.full <- formula
  #TODO(N.Navarro): !!! Need checking formula !!!!
  
  if(is.null(covar)) {
    fm.red <- "y~1"
    ncov=0
  }
  else {
    fm.red <- paste("y~",paste(colnames(covar),collapse="+"),sep="")
    ncov <- ncol(covar)
  }
  
  partial.effect <- matrix(NA,1+ncov+Q$n.qtl*(n.gen-1),length(id.phen))
  for (i in id.phen){
    partial.effect [,i-min(id.phen)+1] <- fitqtl(cross,pheno.col=i,covar=covar,method="hk", qtl=Q,formula=fm.full,dropone=FALSE,get.est=T)$ests$ests
  }
  if (!is.null(pca)){
    partial.effect <- back.rotation(partial.effect,pca)
  } #back to the tangent coordinate space
  rownames(partial.effect) <- names(fitqtl(cross,pheno.col=i,covar=covar,method="hk",qtl=Q,formula=fm.full,get.ests=T)$ests$ests)
  #c("intercept",colnames(covar),rownames(qtls))
  return(partial.effect)
}
#' Effect size for multivariate shape vector
#' 
#' A function to compute the effect size from the result obtained from fitqtl.shape.
#' 
#' Function takes the qtl object, an optionally a shape matrix.
#' 
#' @param qtl A qtl effect obtained from fitqtl.shape
#' @param shape A matrix of shape variables [n.ind n.pheno]
#' @return The function returns the effect size of qtls.
#' @references xxxx
#' @seealso  \code{\link{fitqtl.shape}}
#' @keywords effect size
#' @note If a shape matrix is done, then the function returns a list with the percentage explained of variance of the projection scores. 
#' @author Nicolas Navarro
#' @examples 
#' qtl.effect <- fitqtl.shape(cross, pheno.col=1:2, qtl, covar=NULL, "y~covar")
#' ES <- effect.size()
#' plot.effect(m.shape,target,n.dim,n.land,main="Q1")
#' @export
effectsize.shape <- function(qtl,shape=NULL){
  if(!is.null(shape) && class(shape)!="matrix") stop("shape must be a matrix")
  
  qtl <- qtl[!(rownames(qtl)%in%"Intercept"),]
  ES <- sqrt(rowSums(qtl^2))
  if (!is.null(shape)) {
    SS <- rep(NA,nrow(qtl))
    for (i in 1:nrow(qtl)){
      Reg.proj <- shape %*% qtl[i, ] %*% sqrt(solve(t(qtl[i, ]) %*% qtl[i,]))
      SS[i] <- sum(diag(crossprod(Reg.proj)))	
    }
    perc.SS <- SS/sum(diag(crossprod(scale(shape,scale=F))))*100
    # this is a little different of the % SS predicted
    #predSS <- sum(diag(SSCPfull)) / sum(diag(crossprod(scale(shape,scale=F)))) *100
    names(SS) <- names(perc.SS) <- names(ES)
    
    ES <- list(pD=ES, SS.projScr=SS,'%SS.projScr'=perc.SS)
  }
  return(ES)	
}
#' 2D Lollipop plot of qtl effect
#' 
#' A function to plot scaled qtl effect from the result obtained from fitqtl.shape.
#' 
#' Function takes the mean shape, a target shape, and a scaling factor and plot the shape changes as lollipop graph. The function plot either two- or three-dimensional shapes.
#' 
#' @param mshape The mean shape. An array of [n.landmarks n.dim]
#' @param target An array [n.landmarks n.dim] containing the vector of shape changes
#' @param links An optional matrix defining links between landmarks
#' @seealso  \code{\link{fitqtl.shape}}
#' @keywords plot
#' @author Nicolas Navarro
#' @examples 
#' qtl.effect <- fitqtl.shape(cross, pheno.col=1:2, qtl, covar=NULL, "y~covar")
#' target <- qtl.effect*scaling + m.shape
#' plot.shape(m.shape,target,n.dim,n.land,main="Q1")
#' @export
plot.shape <- function(mshape,target,n.dim,n.land,links,main=NULL,...){
  if (n.dim!=2) par(mfrow=c(2,2))
  for (i in 1:2) 
    for(j in (i+1):n.dim){
      if(n.dim==2&i<2||n.dim!=2){
        plot(mshape[seq(i,n.dim*n.land,by=n.dim)],
             mshape[seq(j,n.dim*n.land,by=n.dim)],
             pch=19,axes=F,xlab="",ylab="",col="black",main=main,...)
        segments(mshape[seq(i,n.dim*n.land,by=n.dim)], 
                 mshape[seq(j,n.dim*n.land,by=n.dim)],
                 x1=target[seq(i,n.dim*n.land,by=n.dim)], 
                 y1=target[seq(j,n.dim*n.land,by=n.dim)], 
                 lwd=2,xlab="",ylab="",col="black")
        text(mshape[seq(i,n.dim*n.land,by=n.dim)], 
             mshape[seq(j,n.dim*n.land,by=n.dim)], 
             labels=c(1:n.land),col="white",cex=.45)
      }}
}
#Back.rotation using non-zero eigenvectors
back.rotation <- function(effect,pca){
  null.eigenvalues <- pca$sdev^2<.Machine$double.eps
  #But people may have selected less variables (eg few 1st PCs)
  if (sum(null.eigenvalues)!=ncol(effect)) null.eigenvalues[(ncol(effect)+1):length(null.eigenvalues)] <- TRUE
  return(effect%*%t(pca$rotation[,!null.eigenvalues]))
}