#' fitqtlShape is a function that allows estimating qtl effect.
#' @title fitqtlShape
#' @author Nicolas Navarro
#' @param cross A cross object.
#' @return The function returns the qtl shape effects, eventually projected back in the tangent space.
#' @references xxxx
#' @keywords shape effect
#' @note If cross is ....
#' @examples
#' eff <- fitqtlShape(cross, pheno.col=1:2, qtl)
#' @export

fitqtlShape <- function(cross, pheno.col, qtl, covar=NULL, formula, 
                        method="hk", model="mvnorm",dropone=FALSE, get.ests=TRUE, run.checks=TRUE, forceXcovar=FALSE, pca=NULL) 
{
    if (!any(class(cross) == "cross")) 
        stop("Cross should have class \"cross\".")
    if (!any(class(cross)%in%c("bc","f2"))) 
        stop("Cross should be an \"f2\" or a \"bc\".")
    if (!any(class(cross) == "shape")) 
        stop("Cross should have class \"shape\".")
    
    #qtl may be q qtl object or a summary from scanone.shape
    if (!any(class(qtl)%in%c("qtl","summary.scanone", "summary.stepwiseqtl"))) 
        stop("qtl argument should have class \"qtl\" or \"summary.scanone\" or \"summary.stepwiseqtl\".")
    if (any(class(qtl)%in%c("summary.scanone", "summary.stepwiseqtl"))) {
        qtls.name <- find.pseudomarker(cross, chr = qtls$chr, pos = qtls$pos, where = "prob")
        Q <- makeqtl(cross, qtls$chr, qtls$pos, qtls.name, what = "prob")
    } else {
        Q <- qtl
    }
    n.gen <- Q$n.gen[1]
    if (any(Q$n.gen != n.gen)) 
        stop("QTLs have different number of alleles")
    
    if(is.character(pheno.col)) {
        if (length(which(colnames(cross$pheno)%in%pheno.col)))
            stop("Phenotype \"", pheno.col, "\" don't match")
        pheno.col <- which(colnames(cross$pheno)%in%pheno.col)
    }
    if(is.null(covar)) {
        fm.red <- "y~1"
        ncov=0
    }
    else {
        fm.red <- paste("y~", paste(colnames(covar), collapse="+"), sep="")
        ncov <- ncol(covar)
    }
    if (missing(formula)) {
        fm.full <- paste(fm.red, paste("Q", 1:Q$n.qtl, sep="", collapse="+"), sep="+")
    } else {
        fm.full <- formula
    }
    #TODO(Nico): !!! Need checking formula !!!!
    
    partial.effect <- matrix(NA, 1+ncov+Q$n.qtl*(n.gen-1), length(pheno.col))
    for (i in pheno.col){
        partial.effect [,i-min(pheno.col)+1] <- fitqtl(cross, pheno.col = i, covar = covar,
                                                       method = "hk", qtl = Q, formula = fm.full,
                                                       dropone = FALSE, get.est = TRUE)$ests$ests
    }
    if (!is.null(pca)){
        partial.effect <- backRotation(partial.effect, pca)
    } #back to the tangent coordinate space
    
    # c("intercept",colnames(covar),rownames(qtls))
    rownames(partial.effect) <- names(fitqtl(cross, pheno.col = i, covar = covar, 
                                             method="hk", qtl = Q, formula = fm.full, 
                                             dropone = FALSE, get.ests = TRUE)$ests$ests)
    
    class(partial.effect) <- c("qtleffect",class(partial.effect))
    return(partial.effect)
}
#' Effect size for multivariate shape vector
#' 
#' Compute the effect size of qtl effect.
#' 
#' Function takes a qtl object, and an optionally a shape matrix.
#' 
#' @param qtl A qtl effect obtained from fitqtlShape
#' @param shape A optional matrix of shape variables [n.ind n.pheno]
#' @param geno A optional matrix of genotype data [n.ind n.qtl]
#' @return The function returns the effect size of qtls.
#' @note If the shape matrix and/or the genotypes are provided, then the function returns a dataframe with the ES standardized by genotype variance and the percentage explained of variance of the projection scores. 
#' @author Nicolas Navarro
#' @examples 
#' ES <- effectsizeShape(qtl, shape, geno = as.matrix(getGeno(cross, qtl)))
#' @export
effectsizeShape <- function(qtl, shape=NULL, geno=NULL, ...){
    if (!is.matrix(qtl)) 
        stop("qtl must be a matrix")
    if (!is.null(shape) && class(shape) != "matrix") 
        stop("shape must be a matrix")

    qtl <- qtl[!(rownames(qtl)%in%"Intercept"), , drop = FALSE]
    ES <- sqrt(rowSums(qtl^2))
    if (!is.null(shape)) {
        SS <- SSprojScres.model <- SSmod <- rep(NA,nrow(qtl))
        if (!is.null(geno)) {
            pD.std <- ES^2*diag(cov(geno))/sum(diag(cov(shape)))
        }
        for (i in 1:nrow(qtl)){
            Reg.proj <- shape %*% qtl[i, ] %*% sqrt(solve(t(qtl[i, ]) %*% qtl[i,]))
            SS[i] <- crossprod(scale(Reg.proj, scale = FALSE))
            if (!is.null(geno)) {
                xx <- as.matrix(geno[, -i])
                mod.red <- lm(Reg.proj ~ xx)
                mod <- lm(Reg.proj ~ geno)
                SSmod[i] <- crossprod(mod.red$residuals) - crossprod(mod$residuals)
                
                plot(Reg.proj ~ geno[,i], xlab = expression(Pr(g[i]==j~"|"~bold(M)[i])),
                     ylab = "proj.Scores", main = rownames(qtl)[i], ...) 
                abline(lm(Reg.proj ~ geno[,i]), ...)
                SSprojScres.resids <- lm(Reg.proj ~ geno[,i])$residuals
                SSprojScres.model[i] <-  SS[i] - crossprod(SSprojScres.resids)
                
            }    
        }
        SST <- sum(diag(crossprod(scale(shape, scale=FALSE))))
        perc.SS <- SSprojScres.model / SST * 100
        # this is a equivalent to the % SS predicted
        # predSS <- sum(diag(SSCPfull)) / sum(diag(crossprod(scale(shape, scale=F)))) *100
        
        # This is the SS of the projection scores given the total amount of variation 
        # How much variance there are in the direction of the qtl effect
        perc.SSprojScres <- SS / SST *100
        # This is the ammount of variance explained in that direction 
        perc.SSprojScres.explained <- SSprojScres.model / SS * 100
        names(SS) <- names(perc.SS) <- names(ES)
        pD <- ES
        perc.SST <- SSmod / SST * 100
        ES <- data.frame(pD, pD.std, perc.SS, perc.SSprojScres, perc.SSprojScres.explained, perc.SST)
    }
    return(ES)    
}  

#' @title 2D Lollipop plot of qtl effect
#' 
#' @description Plots scaled shape changes
#' 
#' @details Function takes the mean shape, a shape change, and a scaling factor and plots the shape changes as lollipop graph. The function plots either two- or three-dimensional shapes.
#' 
#' @param mshape The mean shape. An array of [n.landmarks n.dim]
#' @param effect An array [n.landmarks n.dim] containing the effect (shape changes)
#'  defined as the difference between a target shape and the mean shape
#' @param scaling Scaling factor. The target shape = mean shape + effect * scaling
#' @param links An optional matrix defining links between landmarks
#' @author Nicolas Navarro
#' @examples 
#' x <- matrix(rnorm(10*3*2), nrow=10)
#' x <- asShapeArray(x, n.land = 3, n.dim = 2, byrow = TRUE)
#' mshape <- matrix(c(0,0,1,0,0.5,1), nrow = 3, byrow = TRUE)
#' plot.shapeEffect(mshape, x[, , 1], scaling = 2, main="Triangle 1")
#' @export
plot.shapeEffect <- function(mshape, effect, scaling = 1, 
                             links, labels, col.links, ...) {
    
    #checking
    if (!is.matrix(mshape) && !is.array(mshape)) 
        stop("mshape should be a matrix or an array")
    if (!is.matrix(effect) && !is.array(effect))
        stop("effect should be a matrix or an array")
    if (any(dim(mshape)[1:2] != dim(effect)[1:2]))
        stop("dimensions of inputs didn't agree")
    n.dim <- dim(mshape)[2]
    if (length(dim(mshape)) > 2) mshape <- mshape[, , 1]
    
    # default values:
    cex <- 0.45
    col1 <- "black"
    col2 <- "white"
    if (missing(col.links)) col.links <- "gray"
    pch <- 19
    lwd <- 2
    xlab <- ylab <- ""
    if (missing(labels)) labels <- c(1:n.land)
    # Process suppl args
    argin <- list(...)
    if (length(argin)) {
        if ("col"%in%names(argin)) col1 <- argin$col
        if ("pch"%in%names(argin)) pch <- argin$pch
        if ("lwd"%in%names(argin)) lwd <- argin$lwd
        if ("xlab"%in%names(argin)) xlab<- argin$xlab
        if ("ylab"%in%names(argin)) ylab <- argin$ylab
        if ("labels"%in%names(argin)) labels <- argin$labels
        if ("cex"%in%names(argin)) cex <- argin$cex
    }
    target <- mshape + effect * scaling
    if (n.dim!=2) par(mfrow=c(2,2))
    for (i in 1:2) {
        if(n.dim == 2 & i < 2 || n.dim != 2) {
            plot(c(mshape[, i], target[, i]), c(mshape[, i+1], target[, i+1]), 
                 type = "n", axes = FALSE, xlab = xlab, ylab = ylab, ...)
            if (!missing(links)){
                for (ii in 1:nrow(links)){
                    segments(mshape[links[ii,1], i],mshape[links[ii,1], i+1], 
                             mshape[links[ii,2], i],mshape[links[ii,2], i+1], 
                             lwd = lwd, col = col.links, ...)
                }  
            }
            points(mshape[, i], mshape[, i+1], pch = pch, col = col1)
            segments(mshape[, i], mshape[, i+1],
                     x1 = target[, i], y1 = target[, i+1], lwd = lwd, col = col1)
            text(mshape[, i], mshape[, i+1], 
                 labels = labels, col = col2, cex = cex)
        }
    }
}
#BackRotation using non-zero eigenvectors
backRotation <- function(effect, pca){
    
    if (!is.matrix(effect)) 
        stop("effect must be a matrix [n.effect x n.variables]")
    
    #Check PCA from prcomp
    if (class(pca)=="prcomp"){
        null.eigenvalues <- pca$sdev^2<.Machine$double.eps
        #But people may have selected less variables (eg few 1st PCs)
        if (sum(null.eigenvalues)!=ncol(effect)) 
            null.eigenvalues[(ncol(effect)+1):length(null.eigenvalues)] <- TRUE
        rotation <- pca$rotation[,!null.eigenvalues]
    } else {
        if (!is.matrix(pca))
            stop("pca should be either an object from prcomp() or a matrix of the eigenvectors such as prcomp(X)$rotation ")
        rotation <- pca[,1:ncol(effect)]
    } 
    return(effect%*%t(rotation))
}

#' 
#' @title asShapeArray
#' @author Nicolas Navarro
#' @description Convert a matrix of landmark coordinates into an array [n.land x n.dim x n.obs].
#' @details The function converts a matrix of landmark coordinates into an array. Such arrays are common format used in other morphometrics packages (shape, geomorph, Morpho). Each row of the input matrix contains coordinates for a single specimens.
#' @param shapes A matrix of shape coordinates n.obs rows of [x1 y1 {z1} x2 y2 {z2} ...]
#' @param n.land Number of landmarks
#' @param n.dim Number of dimensions (2 or 3)
#' @param byrow logical (default TRUE: data arranged initially according to x1 y1 z1 x2 y2 z2...). FALSE: data arranged initially according to x1 x2 x3 .... xp y1 .... yp z1 .... zp (used in the Morpho package)
#' @return The function returns a multidimensional array [n.land x n.dim x n.obs]. The third dimension contains the names of each specimens if it was specified as row names in the input matrix.
#' @seealso  \code{\link[geomorph]{arrayspecs}} 
#' @examples
#' x <- matrix(rnorm(5*3*2),nrow=5, byrow=TRUE) # 5 random triangles
#' asShapeArray(x, 3, 2, byrow = TRUE) 
#' @export
asShapeArray <- function(shapes, n.land, n.dim, byrow = TRUE){   
    if (!is.matrix(shapes)) 
        stop("shapes must be a matrix")
    if (prod(dim(shapes)) != nrow(shapes)*n.dim*n.land) 
        stop("dimensions don't match")
    if (!(n.dim%in%c(2,3)))
        stop("Only 2D or 3D allow")
    rowNam <- rownames(shapes)
    if (byrow){
        shapes <- aperm(array(shapes, dim=c(nrow(shapes), n.dim, n.land)),c(3, 2, 1))
    } else {
        shapes <- aperm(array(shapes, dim=c(nrow(shapes), n.land, n.dim)),c(2, 3, 1))
    }
    dimnames(shapes) <- list(NULL,NULL,rowNam)
    return(shapes)  
}
