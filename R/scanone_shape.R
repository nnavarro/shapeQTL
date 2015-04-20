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
                         formula) {
    # TODO(Nico): Check the formula
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
        if (!is.null(addcovar)){
            if (!is.data.frame(addcovar))
                addcovar <- as.data.frame(addcovar)
        }
        if (!is.null(intcovar)) {
            if (!is.data.frame(intcovar)) 
                intcovar <- as.data.frame(intcovar)
            if (is.null(addcovar)) {
                warning("Main effects of covariates added to the model")
                addcovar <- intcovar
            } else {
                isCovarAgree <- names(intcovar)%in%names(addcovar)
                if (any(!isCovarAgree)) {
                    warning("Main effects of covariates added to the model")
                    addcovar <- cbind(addcovar, intcovar[, which(!isCovarAgree)])
                } 
            } 
        }
        if (missing(n.perm)) 
            n.perm <- 0
        if (missing(formula)) {
            if (is.null(addcovar)) {
                formula <- as.formula("pheno ~ 1")
            } else {
                formula <- as.formula(paste("pheno ~", paste(names(addcovar), collapse = "+")))
            }
        } else {
            if (is.character(formula)) 
                formula <- as.formula(formula)
        }
        # Check if missing observations in phenotype
        if (sum(is.na(cross$pheno[,pheno.col])) != 0) {
            warning("Missing observations in phenotypes have been removed")
            isNAp <- (rowSums(is.na(cross$pheno[, pheno.col, drop = FALSE])) != 0)
            cross <- subset(cross, ind = cross$pheno$ID[!isNAp])           
            if (!is.null(addcovar)){
                addcovar <- addcovar[!isNAp, , drop = FALSE]
            }
            if (!is.null(intcovar)){
                intcovar <- intcovar[!isNAp, , drop = FALSE]
            }
        }
        # Check if missing observations in covariates
        if (!is.null(addcovar)){
            if (!is.data.frame(addcovar)) 
                stop("addcovar must be a data.frame")
            isNA <- (rowSums(is.na(addcovar)) != 0)
            if (!is.null(intcovar)){
                if (!is.data.frame(addcovar)) 
                    stop("intcovar must be a data.frame")
                isNA <- isNA + (rowSums(is.na(intcovar) != 0))
            }
            if (sum(isNA) != 0) {
                warning("Missing observations in covariates removed")
                addcovar <- addcovar[!isNA, , drop = FALSE]
                if (!is.null(intcovar)){
                    intcovar <- addcovar[!isNA, , drop = FALSE]
                }
                # Update cross accordingly
                cross <- subset(cross, ind = cross$pheno$ID[!isNA])
            }
        }
        # Checks obs agreement in geno/pheno
        try(nind(cross))
        # Get the multivariate statistics to use - Only one is allowed
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
            genomWideMax[[1]] <- rbind(genomWideMax[[1]], genomWideMax[[j]])
        
        outScan1 <- genomWideMax[[1]]
    } else { 
        if (n.perm > 0 && n.cluster < 2) {
            # If bc result has only 1 lod column, whereas f2 leads to 3 lod columns
            # TODO(Nico): should do more control on cross type?
            genomWideMax <- matrix(0, n.perm, 
                                   ifelse(test = class(cross)[1] == "bc", 
                                          yes = 1, no = 3))
            
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
                genomWideMax[perm, ] <- apply(tmp[,-c(1,2),drop=FALSE], 2, max,
                                              na.rm = TRUE)
            }
            colnames(genomWideMax) <- colnames(tmp)[-c(1,2)]
            outScan1 <- genomWideMax
        } else {
            #---------------------------------------------------
            # 3. One-dimensional scan
            # TODO(Nico): Is pheno must be matrix ?
            outScan1 <- mvGenomScan(cross, pheno = as.matrix(cross$pheno[,pheno.col]),
                                    mod.red = formula, covar = addcovar, 
                                    test = test, chr = chr) # back.qtl = NULL  
        }
    }
    return(outScan1)
}
