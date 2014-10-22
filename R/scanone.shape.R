####################################### 
# shapeQTL mapping experiment with R
#
# Nicolas Navarro - 2013-2014
########################################  
#' Run one-dimensional genomescan 
#' 
#' A function to run QTL mapping on multivariate phenotypes (eg tangent coordinates or PCs)
#' 
#' Function takes a modified cross object with an additional class shape obtained by example with \code{\link{update.cross}}.
#' Additional covariates may be provided as supplementary arguments. Their ordering must match to the new phenotypes.       
#' 
#' @param cross A cross object containing genotypes and some phenotypes
#' @param chr (Optional) Names of chromosomes use for mapping
#' @param pheno.col A numerical vector with column indicies of the phenotypes to map in cross$pheno or a vector containing their names. If missing, All but ID column.
#' @param model The phenotype model. Default and only possible value: \"mvnorm\"
#' @param method Method to use. Default and only value: \"hk\" for Haley-Knott regression
#' @param addcovar A data.frame with the additive covariates
#' @param intcovar A data.frame with the interactive covariates - KEEP IT DEFAULT, not use so far
#' @param use Only \"all.obs\" - KEEP IT DEFAULT, not use so far
#' @param n.perm Number of permutations
#' @param perm.Xsp - KEEP IT DEFAULT, not use so far
#' @param perm.strata - KEEP IT DEFAULT, not use so far
#' @param verbose - KEEP IT DEFAULT, not use so far
#' @param n.cluster Number of clusters for parallelization of permutations
#' @param batchsize - KEEP IT DEFAULT, not use so far
#' @param test Multivariate test statistics: \"Pillai\", \"Hotelling.Lawley\",\"Lik.ratio\". Allows partial matching.
#' @param formula Provides formula for the null model. Optional.
#' @export
#' @seealso  \code{\link{scanone}}
#' @keywords utilitheies
#' @author Nicolas Navarro
#' @return Function returns one dimensional scan simimilar to the R/qtl \code{\link{scanone}} function.
#' @examples 
#' data(fake.bc)
#' scan <- scanoneShape(cross, chr=1, pheno.col=2:3, addcovar=
scanone.shape <- function(cross, chr, pheno.col, model="mvnorm", method="hk", addcovar=NULL, intcovar=NULL,use = "all.obs",
    n.perm, perm.Xsp=FALSE, perm.strata=NULL, verbose=FALSE, batchsize=250, 
    n.cluster = 1,test=c( "Pillai","Hotelling.Lawley","Lik.ratio"),formula){
#covar,mod.red,type="all",back.qtl=NULL,test="Pillai")
	# 1. Checking
	if(!any(class(cross)=="cross")) stop("Cross should have class \"cross\".")
	if(!any(class(cross)%in%c("bc","f2"))) stop("Cross should be an \"f2\" or a \"bc\".")
	if(!any(class(cross)=="shape")) stop("Cross should have class \"shape\".")
	if (batchsize < 1) stop("batchsize must be >=1.")
	if (use!="all.obs") stop("use should be \"all.obs\". Phenoytpes should be complete. Use update.cross to remove missing observations")
	model <- match.arg(model)
	method <- match.arg(method)
	test <- match.arg(test)
  
	if (missing(chr)) chr <- names(cross$geno)
	if (missing(pheno.col)) {
		warning("Missing pheno.col: All columns of cross$pheno but ID taken as phenotypes.")
		pheno.col <- -which(colnames(cross$pheno)%in%c("id","iD","Id","ID"))
	}
	else if(is.character(pheno.col)) {
	  if(length(which(colnames(cross$pheno)%in%pheno.col)))
	    stop("Phenotype \"", pheno.col, "\" don't match")
	  pheno.col <- which(colnames(cross$pheno)%in%pheno.col)
	}
	if (is.null(addcovar) && !is.null(intcovar)){
	    warning("Main effects of covariates added to the model")
	    addcovar <- intcovar
	  }
	if (missing(formula)) {
    if (is.null(addcovar)) formula <- as.formula("pheno ~ 1")
    else formula <- as.formula("pheno ~ covar")
	}
  else if(is.character(formula)) formula <- as.formula(formula)
	# TODO(N.Navarro): Check the formula
  
	# 2. permutation handling
	if(missing(n.perm)) n.perm <- 0
	if(n.perm > 0 && n.cluster > 1) {
		# Some code here are from R/qtl
    	cat(" -Running permutations via a cluster of", n.cluster, "nodes.\n")
    	RNGkind("L'Ecuyer-CMRG")
    	cl <- makeCluster(n.cluster)
    	clusterStopped <- FALSE
    	on.exit(if(!clusterStopped) stopCluster(cl))
    	clusterEvalQ(cl, require(shapeQTL, quietly=TRUE))
    	n.perm <- ceiling(n.perm/n.cluster)

    	genomewide.max <- clusterCall(cl, scanone.shape, cross, chr, 
    				pheno.col, model, method, addcovar, intcovar,
    				use,n.perm, perm.Xsp, 
    				perm.strata, verbose, batchsize, n.cluster = 0,test,formula)				              
    	stopCluster(cl)
    	clusterStopped <- TRUE
    	
    	for(j in 2:length(genomewide.max))
      		genomewide.max[[1]] <- c(genomewide.max[[1]], genomewide.max[[j]])
    	
    	return(genomewide.max[[1]])
  	}
  	else if(n.perm > 0 && n.cluster==0){
		genomewide.max <- rep(0,n.perm)
		n.ind <- nrow(cross$pheno)
		for ( perm in 1:n.perm){
			f <- sample(1:n.ind)
            cross$pheno <- cross$pheno[f, , drop = FALSE]
			if (!is.null(addcovar)) 
                addcovar <- addcovar[f, , drop = FALSE]
            if (!is.null(intcovar)) 
                intcovar <- intcovar[f, , drop = FALSE]
                
      #get only the max -logp values over each chromosome 
			tmp <- scanone.shape(cross, chr, 
    				pheno.col, model, method, addcovar, intcovar,
    				use,n.perm=0, perm.Xsp, 
    				perm.strata, verbose, batchsize, n.cluster = 0, test, formula) 
			genomewide.max[perm] <- apply(tmp[,-c(1,2),drop=FALSE],2, max,na.rm=TRUE)
		}
		return(genomewide.max)
	}
	
	# 3. One-dimensional genome scan
	mv.scan <- multiv.scanone(cross,as.matrix(cross$pheno[,pheno.col]),formula,as.matrix(addcovar),back.qtl=NULL,test)
	return(mv.scan)
}
			
