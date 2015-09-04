#' @title update.cross
#' @description Update phenotypes and covariates in a cross qtl object. New phenotypes (for ex. tangent coordinates) are merged to an existing cross object obtained from \code{\link{read.cross}}
#' @details Function takes a cross object and a dataframe with new phenotypes (tangent coordinates, PC scores,...), matches cross$pheno$ID to the id column of the new phenotypes. Matching id.geno may be provide as an 
#' optional input parameter.
#' Additional covariates may be provided as supplementary arguments but their ordering must match to the new phenotypes.       
#' 
#' @param cross A cross object containing genotypes and some phenotypes
#' @param new.pheno An optional data.frame containing new phenotypes
#' @param phen2keep An optional vector with names of the original phenotypes to keep
#' @param phen2update An optional vector with names of the subset of new 
#' phenotypes to keep
#' @param id.geno An optional vector with the matching id in the genotypes
#' @param na.rm An optional logical if individuals with missing phenotypes 
#' should be removed from the cross
#' @export
#' @seealso  \code{\link{read.cross}}
#' @keywords utilities
#' @author Nicolas Navarro
#' @return Function returns the cross object with updated phenotypes.
#' @examples 
#' data(fake.bc)
#' tgCoords <- data.frame(1:nind(fake.bc), rnorm(nind(fake.bc)), rnorm(nind(fake.bc)))
#' colnames(tgCoords) <- c('Id', 'ProcCoord1','ProcCoord2')
#' fake.bc <- update.cross(fake.bc, new.pheno=tgCoords, phen2keep=c("sex","age"), phen2update=colnames(tgCoords)[grep("ProcCoord",colnames(tgCoords))])

update.cross <- function(cross, new.pheno = NULL, phen2keep = NULL, phen2update = NULL, 
                         id.geno = NULL, na.rm = FALSE, ...) {
    if (!any(class(cross) == "cross")) 
        stop("cross must be a R/qtl cross object")
    #--------------------------------------------------------
    # 1.- update the original cross object
    id <- c("id","Id","ID","iD")
    ID <- NULL
    if (any(id%in%colnames(cross$pheno))) {
        ID <- id[id%in%colnames(cross$pheno)]
        ID <- cross$pheno[, colnames(cross$pheno)%in%ID]
    } else {
        warning("no ID in cross!?: Create an ID = 1 to nind(cross)")
        ID <- 1:nind(cross)
    }	
    if (any(phen2keep%in%id)) 
        phen2keep <- phen2keep[-c(phen2keep%in%id)]
    idx <- match(phen2keep, colnames(cross$pheno), nomatch=0)
    idx <- idx[idx>0]
    if (!length(idx)) {
        cross$pheno <- data.frame(ID = ID)
    } else{  
        cross$pheno <- data.frame(ID = ID, cross$pheno[,idx])
        colnames(cross$pheno) <- c("ID", phen2keep)
    }
    clnms <- colnames(cross$pheno)
    #--------------------------------------------------------
    # 2.- now update the phenotypes	
    if (!is.null(new.pheno)) {
        new.pheno <- data.frame(new.pheno)
        if (is.null(id.geno)) {
            id.geno <- match(cross$pheno[, colnames(cross$pheno)%in%id], 
                             new.pheno[, colnames(new.pheno)%in%id])
        }
        new.pheno <- new.pheno[id.geno, ,drop = FALSE]
        if (!is.null(phen2update)) {
            new.pheno <- new.pheno[, colnames(new.pheno)%in%phen2update, drop = FALSE]
        }
        new.pheno <- as.data.frame(new.pheno)
        clnms <- c(clnms, colnames(new.pheno))
        cross$pheno <- data.frame(cross$pheno, new.pheno)
    }
    #--------------------------------------------------------
    # 3.- process additional covariates
    args <- list(...)
    L <- length(args)
    if (!(!L)) {
        for (i in 1:L) {
            cross$pheno <- data.frame(cross$pheno, args[[i]][id.geno, ])
            cov.names <- paste("Covar",i,sep="")
            clnms <- c(clnms,paste(cov.names, 1:ifelse(is.null(ncol(args[[i]])), 
                                                       1, ncol(args[[i]])), sep="."))
        }
    }
    colnames(cross$pheno) <- clnms
    #--------------------------------------------------------
    # 4.- process missing values in phenotypes
    if (na.rm) {
        isNA <- unique(which(apply(cross$pheno,2,is.na))%%nind(cross))
        if (length(isNA)) {
            isNA[isNA==0] <- nind(cross)
            id <- cross$pheno$ID[-isNA]
            cross <- subset(cross,ind=id)	
        }
    }
    class(cross) <- c(class(cross),"shape")
    return(cross)
}
#' @title getGeno
#' @description Get genotypes in a cross at qtl positions 
#' @details Function takes a cross object and a data.frame with 'chr' and 'pos' and returns genotype probabilities at the closest position in the genetic map.       
#' 
#' @param cross A cross object containing genotypes and some phenotypes
#' @param Q A data.frame containing chr and pos (similar to the output of summary.scanone)
#' @param add.only An optional argument (Default: FALSE). If f2 cross, add.only=TRUE forces the output to be the expected number of B alleles (pure additive model) 
#' @export
#' @keywords utilities
#' @author Nicolas Navarro
#' @return Function returns the genotype with updated phenotypes.
#' @examples 
#' data(fake.bc)
#' fake.bc <- calc.genoprob(fake.bc)
#' Q <- data.frame(chr=1, pos=26)
#' geno <- getGeno(fake.bc, Q)

getGeno <- function(cross, Q, add.only=FALSE){
    if (!any(class(cross)=="cross")) 
        stop("cross must be a R/qtl cross object.")
    if (!any(class(cross)%in%c("bc","f2", "happy"))) 
        stop("this method is implemented only for bc, f2 or happy.") 
    if (!any(class(Q)=="data.frame") || !any(colnames(Q)%in%c("chr","pos"))) 
        stop("Q must be a data.frame with at least a Q$chr and Q$pos columns.")
    if (!("prob"%in%names(cross$geno[[1]]))) 
        stop("You must first run calc.genoprob.")
    if (!("map"%in%names(attributes(cross$geno[[1]]$prob))))
        stop("There is no map associated with the genoprob?")
    
    n.qtl <- nrow(Q)
    n.gen <- rep(2,n.qtl)
    geno <- NULL
    for (q in 1:n.qtl){
        tmp <- subset(cross, chr = Q$chr[q])
        map <- attr(tmp$geno[[1]]$prob, "map")
        q1 <- which.min(abs(map - Q$pos[q]))
        pr <- tmp$geno[[1]]$prob[, q1, , drop = TRUE]
        n.gen[q] <- dim(pr)[2]
        if (class(cross)[1]=="bc" || ((class(cross)[1]=="f2" & !add.only))) {
            pr <- pr[, -dim(pr)[2]]
            n.gen[q] <- n.gen[q] - 1
        } else {
            if(class(cross)[1]=="f2" & add.only) {
                pr <- pr[, 2] + 2*pr[, 3] 
                n.gen[q] <- 1
            }	
        } # we don't remove the last one for happy (handle by QR anyway)
        geno <- cbind(geno,pr)
    }
    gen <- NULL
    for (q in 1:n.qtl) 
        gen <- c(gen, paste(paste("q",q,sep=""), 1:n.gen[q], sep="."))
    colnames(geno) <- gen	
    geno <- as.data.frame(geno)
    return(geno)
}
