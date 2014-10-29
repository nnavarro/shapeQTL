#' fitqtl.shape is a function that allows estimating qtl effect.
#' @title fitqtlShape
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
  if (!any(class(qtl)%in%c("qtl","summary.scanone"))) 
      stop("qtl argument should have class \"qtl\" or \"summary.scanone\".")
  if (any(class(qtl)%in%"summary.scanone")){
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

  
  partial.effect <- matrix(NA, 1+ncov+Q$n.qtl*(n.gen-1), length(id.phen))
  for (i in id.phen){
    partial.effect [,i-min(id.phen)+1] <- fitqtl(cross, pheno.col = i, covar = covar,
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
#' @param qtl A qtl effect obtained from fitqtl.shape
#' @param shape A matrix of shape variables [n.ind n.pheno]
#' @return The function returns the effect size of qtls.
#' @note If a shape matrix is done, then the function returns a list with the percentage explained of variance of the projection scores. 
#' @author Nicolas Navarro
#' @examples 
#' qtl.effect <- fitqtlShape(cross, pheno.col=1:2, qtl, covar=NULL, "y~covar")
#' ES <- effectsizeShape()
#' @export
effectsizeShape <- function(qtl,shape=NULL){
  if (!is.null(shape) && class(shape) != "matrix") stop("shape must be a matrix")
  if (!is.matrix(qtl)) stop("qtl must be a matrix")
  
  qtl <- qtl[!(rownames(qtl)%in%"Intercept"), ]
  ES <- sqrt(rowSums(qtl^2))
  if (!is.null(shape)) {
    SS <- rep(NA,nrow(qtl))
    for (i in 1:nrow(qtl)){
      Reg.proj <- shape %*% qtl[i, ] %*% sqrt(solve(t(qtl[i, ]) %*% qtl[i,]))
      SS[i] <- sum(diag(crossprod(Reg.proj)))	
    }
    perc.SS <- SS / sum(diag(crossprod(scale(shape, scale=FALSE))))*100
    # this is a little different of the % SS predicted
    #predSS <- sum(diag(SSCPfull)) / sum(diag(crossprod(scale(shape, scale=F)))) *100
    names(SS) <- names(perc.SS) <- names(ES)
    
    ES <- list(pD = ES, SS.projScr = SS, '%SS.projScr' = perc.SS)
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
#' x <- asShapeArray(x, n.land = 3, n.dim = 2)
#' mshape <- matrix(c(0,0,1,0,0.5,1), nrow = 3, byrow = TRUE)
#' plot.shapeEffect(mshape, x[, , 1], scaling = 2, main="Triangle 1")
#' @export
plot.shapeEffect <- function(mshape, effect, scaling = 1, 
                             links, labels, col.links, ...) {
                             
    #checking
    if (!is.matrix(mshape) || !is.array(mshape)) 
        stop("mshape should be a matrix or an array")
    if (!is.matrix(effect) || !is.array(effect))
        stop("mshape should be a matrix or an array")
    if (any(dim(mshape) != dim(effect)))
        stop("dimensions of inputs didn't agree")
    n.dim <- dim(mshape)[2]
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
        if (!is.matrix(rotation))
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
#' @return The function returns a multidimensional array [n.land x n.dim x n.obs]. The third dimension contains the names of each specimens if it was specified as row names in the input matrix.
#' @seealso  \code{\link[geomorph]{arrayspecs}} 
#' @examples
#' x <- matrix(rnorm(5*3*2),nrow=5, byrow=TRUE) # 5 random triangles
#' asShapeArray(x, 3, 2) 
#' @export
asShapeArray <- function(shapes, n.land, n.dim){   
    if (!is.matrix(shapes)) 
        stop("shapes must be a matrix")
    if (prod(dim(shapes)) != nrow(shapes)*n.dim*n.land) 
        stop("dimensions don't match")
    if (!(n.dim%in%c(2,3)))
        stop("Only 2D or 3D allow")
    rowNam <- rownames(shapes)
    shapes <- aperm(array(shapes, dim=c(nrow(shapes), n.dim, n.land)),c(3, 2, 1))
    dimnames(shapes) <- list(NULL,NULL,rowNam)
    return(shapes)  
}
