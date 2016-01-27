#' @title fitqtlShape
#' @description Estimates qtl effect 
#' @details Function takes a cross object        
#' 
#' @param cross A cross object containing genotypes and some phenotypes
#' @param pheno.col A numerical vector with column indicies in cross$pheno of 
#' the phenotypes to map or a vector containing their names. 
#' @param qtl An 
#' @param covar An optional data.frame of covariates 
#' @param formula A formula ccc
#' @param method Indicates whether to use multiple imputation or Haley-Knott regression. Only 'hk' is tested so far.
#' @param forceXcovar Not use so far. Keep it default.
#' @param pca An optional prcomp object from \code{\link{prcomp}} function use to get PC scores or just the corresponding eigenvector matrix. Use to return the effect in their original coordinates (PCscores back to tgCoords)
#' @export
#' @keywords utilities
#' @author Nicolas Navarro
#' @return Function returns qtl effects
#' @examples 
#' data(fake.bc)
#' fake.bc <- update.cross(fake.bc, phen2keep=colnames(fake.bc$pheno))
#' fake.bc <- calc.genoprob(fake.bc)
#' qtl <- makeqtl(fake.bc, chr = c(1, 8, 13), pos=c(26, 56, 28), what="prob")
#' eff <- fitqtlShape(fake.bc, pheno.col=grep('pheno',colnames(fake.bc$pheno)), qtl)
#' # Using PC scores
#' pca <- prcomp(fake.bc$pheno[,c('pheno1','pheno2')])
#' PCs <- data.frame(Id=1:nrow(pca$x), pca$x)
#' fake.bc <- update(fake.bc, new.pheno=PCs)
#' eff <- fitqtlShape(fake.bc, pheno.col=grep("PC",colnames(fake.bc$pheno)), qtl, pca=pca)

fitqtlShape <- function(cross, pheno.col, qtl, covar=NULL, formula, 
                        method="hk", forceXcovar=FALSE, pca=NULL) 
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
        qtl.name <- find.pseudomarker(cross, chr = qtl$chr, pos = qtl$pos, where = "prob")
        Q <- makeqtl(cross, qtl$chr, qtl$pos, qtl.name, what = "prob")
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
                                                       method = method, qtl = Q, formula = fm.full,
                                                       dropone = FALSE, get.est = TRUE)$ests$ests
    }
    if (!is.null(pca)){
        partial.effect <- backRotation(partial.effect, pca)
    } #back to the tangent coordinate space
    
    # c("intercept",colnames(covar),rownames(qtls))
    rownames(partial.effect) <- names(fitqtl(cross, pheno.col = i, covar = covar, 
                                             method = method, qtl = Q, formula = fm.full, 
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
#' @param doPlot Logical (default FALSE) to make a plot of the projection scores versus the genotypes
#' @param vcv A list of n.qtl with the [n.pheno n.pheno] dispersion matrix for each qtl
#' @return The function returns the effect size of qtls.
#' @note If the shape matrix and/or the genotypes are provided, then the function returns a dataframe with the ES as the percentage explained of variance of the projection scores. 
#' @author Nicolas Navarro
#' @examples 
#' ES <- effectsizeShape(qtl, shape, geno = as.matrix(getGeno(cross, qtl)))
#' @export
effectsizeShape <- function(qtl, shape=NULL, geno=NULL, doPlot=FALSE, vcv=NULL, ...){
    if (!is.matrix(qtl)) 
        stop("qtl must be a matrix")
    if (!is.null(shape) && class(shape) != "matrix") 
        stop("shape must be a matrix")
    
    qtl <- qtl[!(rownames(qtl)%in%"Intercept"), , drop = FALSE]
    ES <- sqrt(rowSums(qtl^2))
    if (!is.null(shape)) {
        SS <- SSprojScres.model <- SSmod <- perc.SSprojScres.explained <- rep(NA,nrow(qtl))
        names(SS) <- names(ES)
        Z <- scale(shape, scale=FALSE)
        ZZ <- crossprod(Z)
        # total Sum of Squares
        TSS <- sum(diag(ZZ))
        # projected TSS
        SS <- apply(qtl, 1, FUN = projectBetaVCV, SSCP=ZZ)
        if (!is.null(geno)) {
            Reg.proj <- apply(qtl, 1, FUN = projectBetaZ, Z=shape)
            mod <- lm(Reg.proj ~ geno) # estimated multivariate fitted effects
            SSmod <- sapply(1:nrow(qtl), ssmod1, Y = Reg.proj, yhat = mod$fitted.values, g = geno)
        } else if(!is.null(vcv) && is.list(vcv)) {
            if (length(vcv) != nrow(qtl)) stop("length of vcv is different of the number of qtls")
            SSmod <- sapply(1:nrow(qtl), projectBetaVCVlist, qtl, vcv)
        }
        
        # This is the ammount of variance explained in the direction of the qtl effect
        perc.SSprojScres.explained <- SSmod / SS * 100
        # SS of the projection scores given the total amount of variation 
        # How much variance there are in the direction of the qtl effect
        perc.SSprojScres <- SS / TSS *100
        perc.SST <- SSmod / TSS * 100
        ES <- data.frame(pD=ES, perc.SSprojScres, perc.SSprojScres.explained, perc.SST)
        attr(ES, which="TSS") <- TSS
        # Plotting of regression projection scores
        if (doPlot && !is.null(geno)) {
            for (i in 1:nrow(qtl))
                plot(Reg.proj[,i] ~ geno[,i], xlab = expression(Pr(g[i]==j~"|"~bold(M)[i])),
                     ylab = "proj.Scores", main = colnames(X)[i], ...) 
            abline(lm(Reg.proj[,i] ~ geno[,i]), ...)
        }
    }
    return(ES)    
}
# Utility function
projectBetaZ <- function(beta, Z) {
    beta <- matrix(beta, 1, length(beta)) #require to get 1 x k vector
    regProj <- Z %*% t(beta) * 1/sqrt(sum(beta^2))
}
projectBetaVCV <- function(beta, SSCP) {
    beta <- matrix(beta, 1, length(beta))  #require to get 1 x k vector
    SS <- 1/sum(beta^2) * (beta %*% SSCP %*% t(beta))
}
projectBetaVCVlist <- function(i, beta, SSCP) {
    b <- beta[i, , drop=FALSE]
    SS <- 1/sum(b^2) * (b %*% SSCP[[i]] %*% t(b))
}
ssmod1 <- function(i, Y, yhat, g) {
    xx <- as.matrix(g[, -i, drop = FALSE])
    if (ncol(xx) < 1) mod.red <- lm(Y[,i] ~ 1)
    else mod.red <- lm(Y[,i] ~ xx)
    SSmod <- crossprod(yhat[,i] - mod.red$fitted.values)
    return(SSmod)
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
    if (missing(labels)) labels <- c(1:nrow(mshape))
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
# BackRotation using non-zero eigenvectors
# Utility function
backRotation <- function(effect, pca){
    
    if (!is.matrix(effect)) 
        stop("effect must be a matrix [n.effect x n.variables]")
    
    #Check PCA from prcomp
    if (class(pca)=="prcomp"){
        null.eigenvalues <- pca$sdev^2<.Machine$double.eps
        #But people may have selected less variables (eg few 1st PCs)
        if (length(null.eigenvalues)!=ncol(effect)) 
            null.eigenvalues[(ncol(effect)+1):length(null.eigenvalues)] <- TRUE
        rotation <- pca$rotation[,!null.eigenvalues]
    } else {
        if (!is.matrix(pca))
            stop("pca should be either an object from prcomp() or a matrix of the eigenvectors such as prcomp(X)$rotation ")
        rotation <- pca[,1:ncol(effect)]
    } 
    return(effect%*%t(rotation))
}

#' @title asShapeArray
#' @description Convert a matrix of landmark coordinates into an array [n.land x n.dim x n.obs].
#' @details The function converts a matrix of landmark coordinates into an array. Such arrays are common format used in other morphometrics packages (shape, geomorph, Morpho). Each row of the input matrix contains coordinates for a single specimens.
#' @param shapes A matrix of shape coordinates n.obs rows of [x1 y1 {z1} x2 y2 {z2} ...]
#' @param n.land Number of landmarks
#' @param n.dim Number of dimensions (2 or 3)
#' @param byrow logical (default TRUE: data arranged initially according to x1 y1 z1 x2 y2 z2...). FALSE: data arranged initially according to x1 x2 x3 .... xp y1 .... yp z1 .... zp (used in the Morpho package)
#' @export
#' @keywords utilities
#' @author Nicolas Navarro
#' @return The function returns a multidimensional array [n.land x n.dim x n.obs]. The third dimension contains the names of each specimens if it was specified as row names in the input matrix.
#' @seealso  \code{\link[geomorph]{arrayspecs}} 
#' @examples
#' x <- matrix(rnorm(5*3*2),nrow=5, byrow=TRUE) # 5 random triangles
#' asShapeArray(x, 3, 2, byrow = TRUE) 
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
