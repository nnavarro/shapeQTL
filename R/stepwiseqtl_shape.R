#' @title stepwiseqtlShape
#' @description Runs a stepwise multiple QTL model search using penalized negative $log_{10}$ of the p-value
#' @details The function built multiple QTL models using penalized negative $log_{10}$ of the p-value (as a multivariate analogue of the penalized LOD score of the R/qtl package). The search is limited so far to a model without interaction
#' @param cross 
#' @param chr
#' @param pheno.col
#' @param qtl \code{\link[qtl]{summary.scanone}} object or any data.frame with 'chr' 'pos' and 
#'  'lod' values. If missing a one-dimensional scan is run first.
#' @param formula
#' @param max.qtl
#' @param covar
#' @param refine.locations Logical
#' @param additive.only Logical. (Defauft FALSE). Setting it to TRUE may be used 
#'  to force a diplotype model (pure additive) in an intercross (Have to bare in mind
#'  that it is *NOT* as in \code{\link[qtl]{summary.scanone}}). Consistency with lodcolumn
#'  is checked and corrected accordingly.
#' @param scan.pairs Logical. (Default FALSE) #keep it default. Not use so far
#' @param penalties Genomewide threshold 
#' @param distresh numeric (Default 0). Minimal distance, between loci already in the model
#' and a candidate locus, below which the candidate locus is discarded. 
#' @param lodcolumn numeric (1 or 2). With intercross, scanone returns three lod scores (full model, additive only model, and lod for dominance), lodcolumn allows to specify on which one of the full or additive model we want do the QTL selection. 
#' @param keeplodprofile Logical (Default TRUE)
#' @param keeptrace Logical (Default FALSE) 
#' @param verbose Logical (Default FALSE)
#' @param test Multivariate statistics (Default "Pillai". Others are : "Hotelling.Lawley",
#'  "Lik.ratio", "Goodall")
#' @references Broman and Sen 2009. A guide to QTL mapping with R/qtl
#' @export
#' @seealso  \code{\link[qtl]{stepwiseqtl}} 
#' @keywords analysis
#' @author Nicolas Navarro
#' @return The function returns a qtl object.
#' @examples
#' data(fake.bc)
#' fake.bc <- calc.genoprob(fake.bc, step=2.5)
#' fake.bc <- update(fake.bc, phen2keep = names(fake.bc$pheno)) #Just update the class of fake.bc
#' covar <- fake.bc$pheno[, c('sex','age')]
#' out1 <- scanoneShape(fake.bc, pheno.col = 1:2, addcovar = covar,
#' test = "Pillai")
#' Q <- max(out1)
#' outStep <- stepwiseqtlShape(fake.bc, pheno.col = 1:2, qtl = Q, max.qtl = 3, covar = covar, penalties=2.5, verbose = TRUE) 
stepwiseqtlShape <- function(cross, chr, pheno.col = 1, qtl, formula, max.qtl = 10,
                             covar = NULL, 
                             refine.locations = TRUE, 
                             additive.only = FALSE,
                             scan.pairs = FALSE, #keep default
                             penalties,
                             distresh = 0,
                             lodcolumn = 1,
                             keeplodprofile = TRUE,
                             keeptrace = FALSE, 
                             verbose = FALSE,
                             test = "Pillai") {
    # This is the model search version of Broman (p.250)
    pheno <- as.matrix(cross$pheno[, pheno.col])
    if (missing(chr)) 
        chr <- names(cross$geno)
    
    if (missing(penalties))
        stop("you should run first 'scanoneShape' with permutations to get penalties")
    
    if (missing(formula)) {
        formula <- as.formula(paste("pheno ~", paste(names(covar), collapse = " + ")))
    } else if (is.character(formula)) {
        formula <- as.formula(formula)
    }
    
    # Null model (we don't need it except fm.red and except in leave1)
    fm.red <- as.formula(paste("pheno", paste(deparse(formula[-2]), collapse = " + ")))
    mod.red <- lm(fm.red, data = covar)
    SSCPerr.red <- crossprod(mod.red$residuals)    	
    rank.S <- qr(SSCPerr.red)$rank
    
    # Initial model
    if (missing(qtl)) {
        gscan <- scanoneShape(cross, chr = chr, pheno.col = pheno.col, 
                              addcovar = covar, test = test, formula = fm.red)
        qtl <- summary(gscan)
    } else if (class(qtl)[1] != "summary.scanone") {
        if (all(c("chr", "pos", "lod")%in%names(qtl))) {
            qtl <- data.frame(qtl$chr, qtl$pos, qtl$lod)
            names(qtl) <- c("chr", "pos", "lod")
        } else {
            stop("qtl object must have chr, pos as variables")
        }
    }
    #-------------------------------------------------------
    # 0. remove unnessecary lodcolumn and do some checkings on the adequation between
    # additive.only arg and the type of lod scores we want do the selection on
    if (!any(lodcolumn == c(1,2)))
        stop("lodcolumn must be 1 or 2")
    if (ncol(qtl) > 3) {
        qtl <- qtl[, c(1:2, lodcolumn + 2)]
        if (!additive.only & colnames(qtl)[3] == "lod.add") {
            additive.only <- TRUE
            warning("Effects of background QTLs were set to additive only 
                    because selection is done on lod.add")
        } else if (additive.only & colnames(qtl)[3] == "lod") {
            additive.only <- FALSE
            warning("Effects of background QTLs were set to additive+dominance  
                    because selection is done on the full model.")
        }
    }
    colnames(qtl)[3] <- "lod"
    # 1. Select the position with largest LOD score from a single-QTL genome scan
    qtl <- qtl[which.max(qtl$lod), , drop = FALSE]
    n.qtls <- 1
    curplod <- qtl$lod - penalties * n.qtls
    curbestplod <- curplod
    curbest <- qtl
    if (keeptrace) {
        temp <- list(chr = qtl$chr, pos = qtl$pos)
        #attr(temp, "formula") <- deparseQTLformula(formula)
        attr(temp, "pLOD") <- curplod
        class(temp) <- c("compactqtl", "list")
        temp <- list(temp)
        names(temp) <- 0
        thetrace <- c("0", temp)
    }
    #-------------------------------------------------------
    # 2. Forward model search
    iter <- 0
    cat("Forward selection\n")
    while (n.qtls < max.qtl){
        iter <- iter + 1
        back.qtl <- as.matrix(getGeno(cross, Q = qtl, add.only = additive.only))
        n.qtls <- n.qtls + 1
        if (verbose) 
            cat("add qtl:", iter + 1, "\n")
        # 2.1 penalized genome scan
        # mvGenomScan update automatically fm.red + back.qtl if back.qtl is not null
        # but the LOD is not the one of the qtl model vs null (covar only) 
        pLod <- mvGenomScan(cross, pheno, mod.red = fm.red, covar = covar,
                            back.qtl = back.qtl, test = test, chr = chr, 
                            updateRedFormula = FALSE)
        
        # Penalization of Broman and Sen 2009, p.251
        pLod$lod <- pLod$lod - penalties * n.qtls
        # Set to pLod to NA for markers closely linked to QTLs already in the model 
        for (q in 1:(n.qtls-1)) {
            pLod[as.character(pLod$chr) == as.character(qtl$chr[q]) & 
                     abs(pLod$pos - qtl$pos[q]) <= disthresh, lodcolumn + 2] <- NA
        }
        # Get R/qtl max over genome (handle NA)
        tmp <- max(pLod, lodcolumn = lodcolumn)[c(1:2, lodcolumn + 2)]
        colnames(tmp)[3] <- "lod" #ensure colnames consistency after refining positions for rbind()
        # if empty, means no position higher than 0
        # penalty drops LOD scores need updating max.qtl to n.qtls - 1
        if(nrow(tmp) == 0) {
            cat("\n !!!!!!\n
                All pLOD negatives, max.qtl updated to ", n.qtls - 1,
                "\n !!!!!!\n")
            max.qtl <- n.qtls - 1
            n.qtls <- n.qtls - 1
            break
        }
        qtl <- rbind(qtl, tmp) 
        curplod <- tmp$lod 
        if (verbose) {
            cat("new qtl:\n")
            print(qtl[,c(1,2)])
            cat("pLOD = ", curplod, "\n")
        }
        # 2.2 Reorder and refine qtl positions
        # Re-order qtls
        qtl <- qtl[order(match(qtl$chr, chr), qtl$pos), ]
        if (refine.locations) {
            if (verbose) cat("--- Refining positions\n")   
            # 2.2.1 refining qtl positions
            rqtls <- refine.qtl(qtl, cross = cross, fm.red = fm.red, pheno = pheno, 
                                covar = covar, chr = NULL, 
                                max.iter = 10, add.only = additive.only, verbose = FALSE)
            if (any(rqtls$pos != qtl$pos)) { # updated positions
                if (verbose) cat(" ---  Moved a bit\n")
                qtl <- data.frame(rqtls$chr, rqtls$pos, rqtls$lod) #row.names=rqtls$name
                colnames(qtl) <- c("chr", "pos", "lod")
                if(verbose) {
                    cat(" --- New positions\n")
                    print(qtl[,c(1,2)])
                }
                # Updating the penalized LOD score of the model
                geno <- as.matrix(getGeno(cross, Q = qtl, add.only = additive.only))
                tmp.form <- as.formula(paste(paste(deparse(fm.red), collapse = ""), "+ qtl"))
                tmp <- lm.shape.test(geno, pheno = pheno, covar = covar, fm.full = tmp.form,
                                     SSCPerr.red = SSCPerr.red, mod.red.rank = mod.red$rank, 
                                     rank.E = rank.S, test = test)
                
                curplod <- tmp - penalties * n.qtls
                if (verbose) cat("new pLOD = ", curplod, "\n")
            }
            if (verbose) cat("--- Refining done\n")
        }
        
        attr(qtl, "pLOD") <- curplod
        if(curplod > curbestplod) {
            if(verbose){
                cat("** new best ** (pLOD increased by ", 
                    round(curplod - curbestplod, 4), ")\n", sep = "") 
            }
            curbestplod <- curplod
            curbest <- qtl
        }
        if (keeptrace) {
            temp <- list(chr = qtl$chr, pos = qtl$pos)
            #attr(temp, "formula") <- deparseQTLformula(formula)
            attr(temp, "pLOD") <- curplod
            class(temp) <- c("compactqtl", "list")
            temp <- list(temp)
            names(temp) <- iter
            thetrace <- c(thetrace, temp)
        }
        if(n.qtls >= max.qtl) break
    }
    #-------------------------------------------------------
    # 3. Backward model search
    cat("Backward elimination\n")
    while (n.qtls > 1){
        iter <- iter + 1
        if (verbose) cat("remove qtl:",n.qtls,"\n")
        n.qtls <- n.qtls - 1
        tmp <- leave1qtl(cross, qtls = qtl, pheno = pheno, covar = covar, 
                         formula.red = fm.red, mod.red.rank = mod.red$rank, 
                         SSCPerr.red = SSCPerr.red,add.only = additive.only, test = test) 
        # look for best model (the one where the removed qtl has minimal effect)
        tmp$lod <- tmp$lod - penalties * n.qtls
        dropQ <- which.max(tmp$lod) 
        curplod <- max(tmp$lod)
        qtl <- qtl[-dropQ, ] 
        #-------------------------------------------------------    
        # 3.1 Refining qtl positions
        # Re-order qtl (should not make any difference; they should be already ordered) 
        qtl <- qtl[order(match(qtl$chr, chr), qtl$pos), ]
        if (refine.locations & n.qtls > 1) {
            if (verbose) cat("--- Refining positions\n")
            rqtls <- refine.qtl(qtl, cross = cross, fm.red = fm.red, 
                                pheno = pheno, covar = covar, chr = NULL, max.iter = 10)
            if (any(rqtls$pos != qtl$pos)) { # updated positions
                if (verbose) cat(" ---  Moved a bit\n")
                qtl <- data.frame(rqtls$chr, rqtls$pos, rqtls$lod)
                colnames(qtl) <- c("chr", "pos", "lod")
                if(verbose) {
                    cat(" --- New positions\n")
                    print(qtl)
                }
                #-------------------------------------------------------    
                # Updating the penalized LOD score of the model
                geno <- as.matrix(getGeno(cross, Q = qtl, add.only = additive.only))
                tmp.form <- as.formula(paste(paste(deparse(fm.red), collapse =""),"+ qtl"))
                tmp <- lm.shape.test(geno, pheno = pheno, covar = covar, fm.full = tmp.form,
                                     SSCPerr.red = SSCPerr.red, mod.red.rank = mod.red$rank, 
                                     rank.E = rank.S, test = test)
                curplod <- tmp - penalties * n.qtls
            }
        }
        
        attr(qtl, "pLOD") <- curplod
        if(curplod > curbestplod) {
            if(verbose){
                cat("** new best ** (pLOD increased by ", round(curplod - curbestplod, 4),
                    ")\n", sep="")
            }
            curbestplod <- curplod
            curbest <- qtl
        }
        
        if (keeptrace) {
            temp <- list(chr = qtl$chr, pos = qtl$pos)
            #attr(temp, "formula") <- deparseQTLformula(formula)
            attr(temp, "pLOD") <- curplod
            class(temp) <- c("compactqtl", "list")
            temp <- list(temp)
            names(temp) <- iter
            thetrace <- c(thetrace, temp)
        }
    }
    #-------------------------------------------------------
    # 4. drop1qtl the current best model to make sure lod column is SS II partial lod
    qq <- drop1qtl(cross, qtls = curbest, formula.red = fm.red, pheno = pheno, 
                   covar = covar, threshold=0, add.only = additive.only, test = test) 
    curbest$lod <- qq$partial.logp
    rownames(curbest) <- find.pseudomarker(cross, curbest$chr, curbest$pos, where="prob")
    class(curbest) <- c("summary.stepwiseqtl", "data.frame")
    #-------------------------------------------------------
    # 5. Refining positions
    # If so then output is a qtl object instead of a summary
    if (refine.locations) { 
        curbest <- lodprofile.qtl(curbest, cross = cross, fm.red = fm.red, pheno = pheno,
                                  covar = covar, chr = NULL, add.only = additive.only,
                                  test = test)
        if (!keeplodprofile){
            # return only summary
            curbest <- data.frame(chr=tmp$chr,pos=tmp$pos, lod=curbest$lod, row.names = tmp$name)
            class(curbest) <- c("summary.stepwiseqtl","data.frame")
            attr(curbest, which="pLOD") <- curbestplod
        }
    }
    
    if(keeptrace)
        attr(curbest, "trace") <- thetrace
    
    return(curbest)
}
leave1qtl <- function(cross, qtls, pheno, covar, formula.red, mod.red.rank, SSCPerr.red, add.only = FALSE, test = "Pillai", verbose = FALSE){
    # Test a model with one-QTL out versus the null model
    n.qtl <- nrow(qtls)
    n.ind <- nrow(pheno)
    fm.red <- paste(deparse(formula.red[-2]), collapse = "")
    #get genotypes for all qtls
    geno <- getGeno(cross, Q = qtls, add.only = FALSE)
    #fit the full model
    gen <- colnames(geno)
    tmp.dat <- cbind(covar, geno)
    #fit the full minus 1 QTL 
    lod <- rep(0, n.qtl)
    for (q in 1:n.qtl){
        idx <- grep(paste("q", q, ".", sep = ""), gen, fixed=TRUE)
        fm.new <- paste(fm.red, paste(gen[-idx], collapse = "+"), sep="+")
        tmp.x <- model.matrix(as.formula(fm.new), data = tmp.dat)
        mod.full <- .Call(C_CdqrlsShapeQTL, x = tmp.x, y = pheno, tol = 1e-07, chk = FALSE)
        
        dfeff <- mod.full$rank - mod.red.rank
        dferr <- n.ind - mod.full$rank	
        SSCPerr.full <- crossprod(mod.full$residuals)
        rank.E <- qr(SSCPerr.full)$rank
        #partial F-test: Full model vs Reduced model
        SSCPfull <- SSCPerr.red - SSCPerr.full
        if (pmatch(test,"Pillai",nomatch=0)) lod[q] <- Pillai.test(SSCPfull,SSCPerr.full,dfeff,dferr,rank.E)
        if (pmatch(test,"Lik.ratio",nomatch=0)) lod[q] <- LikelihoodRatio.test(SSCPerr.full,SSCPerr.red,dfeff,dferr,rank.E)
        if (pmatch(test,"Hotelling.Lawley",nomatch=0)) lod[q] <- Hotelling.test(SSCPfull,SSCPerr.full,dfeff,dferr,rank.E)
        if (pmatch(test,"GoodallF",nomatch=0)) lod[q] <- goodallF.test(diag(SSCPerr.full), diag(SSCPerr.red), dferr, n.ind - mod.red.rank, rank.E)
    }
    
    qtls <- data.frame(qtls$chr, qtls$pos, lod)
    colnames(qtls) <- c("chr", "pos", "lod") 
    
    return(qtls)
}
drop1qtl <- function(cross,qtls,formula.red,pheno,covar,threshold,add.only=FALSE,test="Pillai",verbose=FALSE){
    cl.qtl <- class(qtls)
    n.qtl <- nrow(qtls)
    if (n.qtl < 2) {
        warning("Less 1 QTL in your model! drop analysis not done")
        partial.logp <- qtls$lod
        names(partial.logp) <- rownames(qtls)
    } else {
        n.ind <- nrow(pheno)
        fm.red <- paste(deparse(formula.red[-2]), collapse = "")
        # get genotypes for all qtls
        geno <- getGeno(cross, Q = qtls, add.only = add.only)
        # fit the full model
        gen <- colnames(geno)
        form <- as.formula(paste(fm.red, paste(gen, collapse = "+"), sep="+"))
        x <- model.matrix(form, data = cbind(covar,geno))
        mod.full <- .Call(C_CdqrlsShapeQTL, x = x, y = pheno, tol = 1e-07, chk = FALSE)
        SSCPerr.full <- crossprod(mod.full$residuals)
        rank.E <- qr(SSCPerr.full)$rank
    
        # fit the full minus 1 QTL to get SSII p-value
        partial.logp <- rep(0, n.qtl)
        for (q in 1:n.qtl){
            idx <- grep(paste("q",q,".",sep=""), gen, fixed=TRUE)
            tmp.x <- model.matrix(as.formula(paste(fm.red, paste(gen[-idx], collapse="+"), sep="+")), data=cbind(covar,geno))
            mod.red <- .Call(C_CdqrlsShapeQTL, x = tmp.x, y = pheno, tol = 1e-07, chk = FALSE)
            dfeff <- mod.full$rank - mod.red$rank
            dferr <- n.ind - mod.full$rank	
            SSCPerr.red <- crossprod(mod.red$residuals)	
        #partial F-test: Full model vs Reduced model
            SSCPfull <- SSCPerr.red - SSCPerr.full
            if (pmatch(test,"Pillai",nomatch=0)) partial.logp[q] <- Pillai.test(SSCPfull,SSCPerr.full,dfeff,dferr,rank.E)
            if (pmatch(test,"Lik.ratio",nomatch=0)) partial.logp[q] <- LikelihoodRatio.test(SSCPerr.full,SSCPerr.red,dfeff,dferr,rank.E)
            if (pmatch(test,"Hotelling.Lawley",nomatch=0)) partial.logp[q] <- Hotelling.test(SSCPfull,SSCPerr.full,dfeff,dferr,rank.E)
            if (pmatch(test,"GoodallF",nomatch=0)) lod[q] <- goodallF.test(diag(SSCPerr.full), diag(SSCPerr.red), dferr, n.ind - mod.red.rank, rank.E)
        }
    }
    names(partial.logp) <- rownames(qtls)
    if (any(partial.logp<threshold)) {
        if (verbose) {
            print("removed qtl")
            print(qtls[partial.logp < threshold, ])
        }
        qtls <- qtls[partial.logp > threshold, ]
        cat("some qtls have been drop!\n updating...\n")
        qtls <- drop1qtl(cross, qtls, formula.red, pheno,covar,threshold,add.only)
    }
    else qtls <- cbind(qtls, partial.logp)
    class(qtls) <- c("summary.drop1qtl","data.frame")
    return(qtls)
}

lodprofile.qtl<- function (qtls, cross, fm.red, pheno, covar, chr = NULL, add.only = FALSE,
                           test = "Pillai", verbose = FALSE) 
{
    if(any(class(qtls)%in%"summary.stepwiseqtl")){
        rownames(qtls) <- qtls$name
        qtls <- qtls[,c("chr","pos","lod")] #suppress name and n.gen to match summary.scanone
    }
    else if(!any(class(qtls)%in%"summary.scanone")) 
        stop("qtls must be either summary.scanone object or summary.stepwiseqtl object")
    
    fm.red <- as.formula(paste("pheno", paste(deparse(fm.red[-2]), collapse = "")))
    mod.red <- lm(fm.red, data = covar)
    if (!is.null(chr)) {
        if (length(chr) > 1) {
            chr <- chr[1]
            cat("Analysis of chr", chr, "only\n")
        }
        back.qtl <- qtls[qtls$chr != chr, ]
        qtls <- qtls[qtls$chr == chr, ]
        if (!length(back.qtl)) 
            break
        if (!length(qtls)) 
            stop("no qtl for this chromosome")
        back.geno <- getGeno(cross, Q = back.qtl, add.only)
        covar <- cbind(covar, back.geno)
        # Update the formula of the reduced model
        fm.red <- as.formula(paste(paste(deparse(fm.red), collapse = ""), 
                                   paste(names(back.geno), collapse = "+"), sep = "+"))
    }
    n.qtls <- nrow(qtls)
    lastOut <- NULL
    newpos <- qtls 
    if (n.qtls == 1)  {
        tmp.cross <- subset(cross, chr = qtls$chr)
        tmp <- mvGenomScan(tmp.cross, pheno = pheno, mod.red = fm.red, 
                           covar = covar, back.qtl = NULL, test = test)
        if (ncol(tmp) > 3) 
            tmp <- tmp[, 1:3]
        tmp[is.na(tmp[, 3]), 3] <- 0
        lastOut[[1]] <- tmp
        newpos[1, ] <- max(tmp)
    } else {    
        for (q in 1:n.qtls) {
            geno <- as.matrix(getGeno(cross, Q = qtls[-q, ], add.only))
            tmp.cross <- subset(cross, chr = qtls[q, "chr"])
            tmp <- mvGenomScan(tmp.cross, pheno = pheno, mod.red = fm.red, 
                           covar = covar, back.qtl = geno, test = test)
            if (ncol(tmp) > 3) 
            tmp <- tmp[, 1:3]
            tmp[is.na(tmp[, 3]), 3] <- 0
            lastOut[[q]] <- tmp
            newpos[q, ] <- max(tmp)
        }
    }
    qtls <- newpos
    qtls.name <- find.pseudomarker(cross, qtls$chr, qtls$pos, where = "prob")
    QTL <- makeqtl(cross, qtls$chr, qtls$pos, qtls.name, what = "prob")
    QTL$lod <- qtls$lod
    names(lastOut) <- QTL$name
    attr(QTL, "lodprofile") <- lastOut
    return(QTL)
}
refine.qtl <- function(qtls, cross, fm.red, pheno, covar, 
                       chr = NULL, 
                       max.iter = 10, 
                       add.only = FALSE,
                       verbose = FALSE, 
                       test = "Pillai") {
    
    fm.red <- as.formula(paste("pheno", paste(deparse(fm.red[-2]), collapse = "")))
    mod.red <- lm(fm.red, data = covar)
    # In the process, we don't update mod.red because it corresponds initially to pheno ~ covar 
    # and the prg if back.qtl update it to pheno ~ covar + back.qtl
    
    run.model <- TRUE
    iter <- 0
    if (!is.null(chr)) {
        if (length(chr) > 1) {
            chr <- chr[1]
            cat("Analysis of chr", chr, "only\n")
        }
        back.qtl <- qtls[qtls$chr != chr,]
        qtls <- qtls[qtls$chr == chr, ]
        if (!length(back.qtl)) break
        if (!length(qtls)) stop("no qtl for this chromosome")
        back.geno <- getGeno(cross, Q = back.qtl, add.only)
        covar <- cbind(covar, back.geno)
        # Update the formula of the reduced model
        fm.red <- as.formula(paste(paste(deparse(fm.red), collapse = ""), 
                                   paste(names(back.geno), collapse = "+"), sep = "+"))
    }
    
    qtls.name <- find.pseudomarker(cross, qtls$chr, qtls$pos, where="prob") 
    QTL <- makeqtl(cross, qtls$chr, qtls$pos, qtls.name, what="prob")
    n.qtls <- QTL$n.qtl
    
    if (is.null(chr) & n.qtls == 1) {
        warning("Only one qtl without background qtl: No refining done")
        return(QTL)
    }
    
    map <- attr(QTL, "map")
    if (is.null(map)) 
        stop("Input qtl object should contain the genetic map.")
    mind <- min(sapply(map, function(a) {
        if (is.matrix(a)) a <- a[1, ]
        min(diff(a))
    }))/2
    if (mind <= 0) mind <- 1e-06
    target.qtl <- NULL
    while (run.model){
        iter <- iter + 1
        old.qtls <- qtls
        # Randomly pick a qtl. But make sure the 1st was not the first before
        o <- sample(n.qtls)
        if(!is.null(target.qtl))
            while(o[1] == target.qtl[1]) o <- sample(n.qtls)
        target.qtl <- o
        
        outStep <- NULL
        if(verbose) cat("------ scan ...")
        for (q in 1:n.qtls) {
            geno <- as.matrix(getGeno(cross, Q = qtls[-target.qtl[q], ], add.only))
            tmp.cross <- subset(cross, chr = qtls[target.qtl[q], 'chr'])
            tmp <- mvGenomScan(tmp.cross, pheno, mod.red = fm.red, covar = covar,
                               back.qtl = geno, test = test)
            # So far, selection only on full model (not just on dominance or additive only) 
            if (ncol(tmp) > 3) { 
                tmp <- tmp[, 1:3, drop = FALSE]
            }
            tmp[is.na(tmp[, 3]), 3] <- 0
            mx.tmp <- max(tmp) #get R/qtl max over the chromosome
            if (verbose) {
                print("back qtl")
                print(qtls[, 1:3])
                print("new qtl")
                print(mx.tmp)
            }
            qtls[target.qtl[q], ] <- mx.tmp
            qtls[target.qtl[q], 'chr'] <- as.character(mx.tmp$chr) #this is because of chr is a factor
            outStep[[q]] <- tmp
        }
        if(verbose) cat(" done\n")
        #check if positions have changed.
        if (max(abs(old.qtls[target.qtl, 'pos']-qtls[target.qtl, 'pos'])) < mind){
            converged <- TRUE
            break
        }
        if (iter>=max.iter) break
    }
    idx <- order(as.numeric(qtls$chr) * max(qtls$pos) + qtls$pos)
    qtls <- qtls[idx, ]
    lastOut <- NULL
    newpos <- qtls
    for (q in 1:n.qtls){
        geno <- as.matrix(getGeno(cross, qtls[-q, ], add.only))
        tmp.cross <- subset(cross, chr = qtls[q, 'chr'])
        tmp <- mvGenomScan(tmp.cross, pheno, mod.red = fm.red, covar = covar,
                           back.qtl = geno, test = test)
        if (ncol(tmp) > 3) { 
            tmp <- tmp[, 1:3]
        }  
        tmp[is.na(tmp[, 3]), 3] <- 0
        lastOut[[q]] <- tmp
        mx.tmp <- max(tmp)
        qtls[q, ] <- mx.tmp
        qtls[q, 'chr'] <- as.character(mx.tmp$chr) #!! because chr factor
    }
    qtls.name <- find.pseudomarker(cross, qtls$chr, qtls$pos, where = "prob")	
    QTL <- makeqtl(cross, qtls$chr, qtls$pos, qtls.name, what = "prob")
    # add the current partial lod scores
    QTL$lod <- qtls$lod
    # add the partial lod score profiles
    names(lastOut) <- QTL$name
    attr(QTL, "lodprofile") <- lastOut 
    
    return(QTL)
}