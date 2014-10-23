####################################### 
# shapeQTL mapping experiment with R
#
# Nicolas Navarro - 2013-2014
########################################  
#' Run one-dimensional genomescan 
#' 
#' A function to run QTL mapping on multivariate phenotypes 
#' 
#' Function takes a modified cross object with an additional class shape 
#' obtained from \code{\link{update.cross}}.
#' Additional covariates may be provided as supplementary arguments. Ordering of 
#' observations in those additional covariates must match the new phenotypes.       
#' 
#' @param cross A cross object containing genotypes and some phenotypes
#' @param chr (Optional) Names of chromosomes use for mapping
#' @param pheno.col A numerical vector with column indicies in cross$pheno of 
#' the phenotypes to map or a vector containing their names. 
#' If missing, All but \"ID\" columns will be chosen as phenotypes.
#' @param addcovar A dataframe with the additive covariates
#' @param intcovar A dataframe with the interactive covariates - KEEP IT DEFAULT
#' This argument and the followings except if state otherwise are not use so far
#' They are there for direct
#' @param n.perm Number of permutations
#' @param perm.Xsp - KEEP IT DEFAULT, not use so far
#' @param perm.strata - KEEP IT DEFAULT, not use so far
#' @param n.cluster Number of clusters for parallelization of permutations
#' @param test Multivariate test statistics: \"Pillai\", \"Hotelling.Lawley\",
#' \"Lik.ratio\". Allows partial matching.
#' @param formula Provides formula for the null model. Optional.
#' @export
#' @seealso  \code{\link{scanone}}
#' @keywords qtl, shape
#' @author Nicolas Navarro
#' @return Function returns one dimensional scan similar to the R/qtl \code{\link{scanone}} function.
#' @examples 
#' data(fake.bc)
#' out <- scanoneShape(cross, chr = 1, pheno.col = 1:2, addcovar = getsex(cross)$sex,
#' test = "Pillai") 
scanoneShape <- function(cross, chr, pheno.col, 
                          addcovar      = NULL, 
                          intcovar      = NULL, 
                          n.perm, 
                          perm.Xsp      = FALSE, 
                          perm.strata   = NULL, 
                          n.cluster     = 1, 
                          test          = c("Pillai","Hotelling.Lawley","Lik.ratio"),
                          formula){
    # TODO(Nico): Check the formula
    # TODO(Nico): Handle pheno ~ sex + logCS + ...
    # TODO(Nico): check what happens when several test name...
    #---------------------------------------------------
	# 1. Error checking (doesn't happen if we are in the cluster loop)
	if (n.cluster > -1) {
	    if (!any(class(cross)=="cross")) 
	        stop("Cross should have class \"cross\".")
	    if (!any(class(cross)%in%c("bc","f2"))) 
	        stop("Cross should be an \"f2\" or a \"bc\".")
	    if (!any(class(cross)=="shape")) 
	        stop("Cross should have class \"shape\".")
	    if (missing(chr)) 
	        chr <- names(cross$geno)
	    if ("X"%in%chr) 
            warning("in case of the X chromosome, mapping may be wrong")
	    if (missing(pheno.col)) {
	        warning("Missing pheno.col: All columns but 'ID' taken as phenotypes.")
	        pheno.col <- -which(colnames(cross$pheno)%in%c("id","iD","Id","ID"))
	    } else {
	        if (is.character(pheno.col)) {
	            if (length(which(colnames(cross$pheno)%in%pheno.col)))
	                stop("Phenotype \"", pheno.col, "\" don't match")
	            pheno.col <- which(colnames(cross$pheno)%in%pheno.col)
	        }
	    }
	    if (is.null(addcovar) && !is.null(intcovar)){
	        warning("Main effects of covariates added to the model")
	        addcovar <- intcovar
	    }
	    if (missing(n.perm)) 
	        n.perm <- 0
	    if (missing(formula)) {
	        if (is.null(addcovar)) {
	            formula <- as.formula("pheno ~ 1")
	        } else {
	            formula <- as.formula("pheno ~ covar")
	        }
	    } else {
	        if(is.character(formula)) 
	            formula <- as.formula(formula)
	    }
	    if (length(apply(cross$pheno[,pheno.col],c(1,2),is.na))) {
	        warning("Missing observations in phenotypes have been removed")
	        cross <- update.cross(cross, new.pheno = cross$pheno, na.rm = TRUE)
	    }
	    # Checks obs agreement in geno/pheno
	    try(nind(cross))
	    # Get the multivariate statistics to use
        test <- match.arg(test)
	}
    
    #---------------------------------------------------
    # 2. permutation handling
	if (n.perm > 0 && n.cluster > 1) {
		# Some code here are from R/qtl
    	cat(" -Running permutations via a cluster of", n.cluster, "nodes.\n")
    	RNGkind("L'Ecuyer-CMRG")
    	cl <- makeCluster(n.cluster)
    	clusterStopped <- FALSE
    	on.exit(if(!clusterStopped) stopCluster(cl))
    	clusterEvalQ(cl, require(shapeQTL, quietly=TRUE))
    	n.perm <- ceiling(n.perm/n.cluster)

    	genomWideMax <- clusterCall(cl, scanoneShape, 
                                cross, chr, pheno.col, addcovar, intcovar,
    				            n.perm, perm.Xsp, perm.strata, n.cluster = -1,
                                test, formula)				              
    	stopCluster(cl)
    	clusterStopped <- TRUE
    	
    	for (j in 2:length(genomWideMax))
    	    genomWideMax[[1]] <- c(genomWideMax[[1]], genomWideMax[[j]])
    	
    	out <- genomWideMax[[1]]
  	} else { 
          if (n.perm > 0 && n.cluster < 2) {
              genomWideMax <- rep(0, n.perm)
              n.ind <- nrow(cross$pheno)
              for (perm in 1:n.perm) {
                  f <- sample(1:n.ind)
                  cross$pheno <- cross$pheno[f, , drop = FALSE]
                  if (!is.null(addcovar)) 
                      addcovar <- addcovar[f, , drop = FALSE]
                  if (!is.null(intcovar)) 
                      intcovar <- intcovar[f, , drop = FALSE]
                  # Get only the max -logp values over each chromosome 
                  tmp <- scanoneShape(cross, chr, 
    				    pheno.col, addcovar, intcovar,
    				    n.perm=0, perm.Xsp, perm.strata, n.cluster = -1, 
                        test, formula) 
                  genomWideMax[perm] <- apply(tmp[,-c(1,2),drop=FALSE],2, max, na.rm=TRUE)
		        }
              out <- genomWideMax
          } else {
              #---------------------------------------------------
              # 3. One-dimensional scan
              mvScan1 <- mvGenomScan(cross, as.matrix(cross$pheno[,pheno.col]),
                                     formula, as.matrix(addcovar), 
                                     back.qtl = NULL, test, chr)
              out <- mvScan1  
          }
  	}
    return(out)
}
			
