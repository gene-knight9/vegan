`anova_lm` <- function (object, ...) 
{
    if (length(list(object, ...)) > 1L) 
        return(stats:::anova.lmlist(object, ...))
    if (!inherits(object, "lm")) 
        warning("calling anova.lm(<fake-lm-object>) ...")
    w <- object$weights
    ssr <- sum(if (is.null(w)) object$residuals^2 else w * object$residuals^2)
    mss <- sum(if (is.null(w)) object$fitted.values^2 else w * 
        object$fitted.values^2)
    if (ssr < 1e-10 * mss) 
        warning("ANOVA F-tests on an essentially perfect fit are unreliable")
    dfr <- df.residual(object)
    p <- object$rank
    if (p > 0L) {
        p1 <- 1L:p
        comp <- object$effects[p1]
        asgn <- object$assign[stats:::qr.lm(object)$pivot][p1]
        nmeffects <- c("(Intercept)", attr(object$terms, "term.labels"))
        tlabels <- nmeffects[1 + unique(asgn)]
        ss <- c(vapply(split(comp^2, asgn), sum, 1), ssr)
        df <- c(lengths(split(asgn, asgn)), dfr)
    }
    else {
        ss <- ssr
        df <- dfr
        tlabels <- character()
    }
    ms <- ss/df
    f <- ms/(ssr/dfr)
    P <- pf(f, df, dfr, lower.tail = FALSE)
    table <- data.frame(df, ss, ms, f, P)
    table$P <- format(table$P, scientific = TRUE, digits = 8)
    table[length(P), 4:5] <- NA
    dimnames(table) <- list(c(tlabels, "Residuals"), c("Df", 
        "Sum Sq", "Mean Sq", "F value", "Pr(>F)"))
    if (attr(object$terms, "intercept")) 
        table <- table[-1, ]
    return(table)
    # structure(table, heading = c("Analysis of Variance Table\n", 
    #     paste("Response:", deparse(formula(object)[[2L]]))), 
    #     class = c("anova", "data.frame"))
}

`adonis2` <-
    function(formula, data, permutations = 999, method = "bray",
             sqrt.dist = FALSE, add = FALSE, by = NULL,
             parallel = getOption("mc.cores"), na.action = na.fail,
             strata = NULL, ...)
{
    ## handle missing data
    if (missing(data))
        data <- model.frame(delete.response(terms(formula)),
                            na.action = na.action)
    ## we accept only by = "terms", "margin", "onedf" or NULL
    if (!is.null(by))
        by <- match.arg(by, c("terms", "margin", "onedf"))
    ## evaluate lhs
    YVAR <- formula[[2]]
    lhs <- eval(YVAR, parent.frame(), environment(formula))
    environment(formula) <- environment()
    ## Take care that input lhs are dissimilarities
    if ((is.matrix(lhs) || is.data.frame(lhs)) &&
        isSymmetric(unname(as.matrix(lhs))))
        lhs <- as.dist(lhs)
    if (!inherits(lhs, "dist"))
        lhs <- vegdist(as.matrix(lhs), method=method, ...)
    ## adjust distances if requested
    if (sqrt.dist)
        lhs <- sqrt(lhs)
    if (is.logical(add) && add)
        add <- "lingoes"
    if (is.character(add)) {
        add <- match.arg(add, c("lingoes", "cailliez"))
        if (add == "lingoes") {
            ac <- addLingoes(as.matrix(lhs))
            lhs <- sqrt(lhs^2 + 2 * ac)
        }
        else if (add == "cailliez") {
            ac <- addCailliez(as.matrix(lhs))
            lhs <- lhs + ac
        }
    }
    ## adonis0 & anova.cca should see only dissimilarities (lhs)
    if (!missing(data)) # expand and check terms
        formula <- terms(formula, data=data)
    if (is.null(attr(data, "terms"))) # not yet a model.frame?
        data <- model.frame(delete.response(terms(formula)), data,
                            na.action = na.action)
    formula <- update(formula, lhs ~ .)
    sol <- adonis0(formula, data = data, method = method)
    ## handle permutations
    perm <- getPermuteMatrix(permutations, NROW(data), strata = strata)
    out <- anova_lm(sol, permutations = perm, by = by,
                 parallel = parallel)
    ## attributes will be lost when adding a new column
    att <- attributes(out)
    ## add traditional adonis output on R2
    out <- rbind(out, "Total" = c(nobs(sol)-1, sol$tot.chi, NA, NA))
    out <- cbind(out[,1:2], "R2" = out[,2]/sol$tot.chi, out[,3:4])
    ## Fix output header to show the adonis2() call instead of adonis0()
    att$heading[2] <- deparse(match.call(), width.cutoff = 500L)
    att$names <- names(out)
    att$row.names <- rownames(out)
    attributes(out) <- att
    out
}

`adonis0` <-
    function(formula, data=NULL, method="bray")
{
    ## First we collect info for the uppermost level of the analysed
    ## object
    Trms <- terms(data)
    sol <- list(call = match.call(),
                method = "adonis",
                terms = Trms,
                terminfo = list(terms = Trms))
    sol$call$formula <- formula(Trms)
    TOL <- 1e-7
    lhs <- formula[[2]]
    lhs <- eval(lhs, environment(formula)) # to force evaluation
    formula[[2]] <- NULL                # to remove the lhs
    rhs <- model.matrix(formula, data) # and finally the model.matrix
    assign <- attr(rhs, "assign") ## assign attribute
    sol$terminfo$assign <- assign[assign > 0]
    rhs <- rhs[,-1, drop=FALSE] # remove the (Intercept) to get rank right
    rhs <- scale(rhs, scale = FALSE, center = TRUE) # center
    qrhs <- qr(rhs)
    ## input lhs should always be dissimilarities
    if (!inherits(lhs, "dist"))
        stop("internal error: contact developers")
    if (any(lhs < -TOL))
        stop("dissimilarities must be non-negative")
    ## if there was an na.action for rhs, we must remove the same rows
    ## and columns from the lhs (initDBRDA later will work similarly
    ## for distances and matrices of distances).
    if (!is.null(nas <- na.action(data))) {
        lhs <- as.matrix(lhs)[-nas,-nas, drop=FALSE]
        n <- nrow(lhs)
    } else
        n <- attr(lhs, "Size")
    ## G is -dmat/2 centred
    G <- initDBRDA(lhs)
    ## preliminaries are over: start working
    Gfit <- qr.fitted(qrhs, G)
    Gres <- qr.resid(qrhs, G)
    ## collect data for the fit
    if(!is.null(qrhs$rank) && qrhs$rank > 0)
        CCA <- list(rank = qrhs$rank,
                    qrank = qrhs$rank,
                    tot.chi = sum(diag(Gfit)),
                    QR = qrhs)
    else
        CCA <- NULL # empty model
    ## collect data for the residuals
    CA <- list(rank = n - max(qrhs$rank, 0) - 1,
               u = matrix(0, nrow=n),
               tot.chi = sum(diag(Gres)))
    ## all together
    sol$tot.chi <- sum(diag(G))
    sol$adjust <- 1
    sol$Ybar <- G
    sol$CCA <- CCA
    sol$CA <- CA
    class(sol) <- c("adonis2", "dbrda", "rda", "cca")
    sol
}
