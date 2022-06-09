#' @title Perform PaIRKAT on the output from the GatherNetworks function
#' @name PaIRKAT
#' @description
#'
#' Pathway Integrated Regression-based Kernel Association Test (PaIRKAT) is a
#' model framework for assessing statistical relationships between networks
#' and some outcome of interest while adjusting for
#' potential confounders and covariates.
#'
#' Use of PaIRKAT is motivated by the analysis of networks of metabolites
#' from a metabolomics assay and the relationship of those networks with a
#' phenotype or clinical outcome of interest, though the method can be
#' generalized to other domains.
#'
#' @param formula.H0 The null model in the "formula" format used
#' in \code{\link[stats]{lm}} and \code{\link[stats]{glm}} functions
#' @param networks networks object obtained
#' with \code{\link[pairkat]{GatherNetworks}}
#' @param tau A parameter to control the amount of smoothing,
#' analagous to a bandwidth parameter in kernel smoothing. We found 1 often
#' gave reasonable results, as over-smoothing can lead to
#' inflated Type I errors.
#'
#' @details The PaIRKAT method is to update the feature matrix, \eqn{Z},
#' with the regularized normalized Laplacian, \eqn{L_R}, before performing
#' the kernel association test. \eqn{L_R} is calculated using a "linear"
#' regularization, \deqn{L_R = (I +\tau L)^-1,} where \eqn{I} is the identity
#' matrix, \eqn{\tau} is a regularization parameter that controls the amount
#' of smoothing, and \eqn{L} is the graph's normalized Laplacian. The updated
#' feature matrix, \eqn{Z*L_R} is matrix used for the kernel association
#' test. \cr
#' The linear regularization and Gaussian kernel is used for all tests.
#' See Carpenter 2021 for complete details on PaIRKAT and Smola 2003
#' for information about graph regularization
#'
#' @references
#' Carpenter CM, Zhang W, Gillenwater L, Severn C, Ghosh T, Bowler R, et al.
#' PaIRKAT: A pathway integrated regression-based kernel association test
#' with applications to metabolomics and COPD phenotypes.
#' bioRxiv. 2021 Apr 26;2021.04.23.440821.
#'
#' Smola AJ, Kondor R. Kernels and Regularization on Graphs.
#' In: Schölkopf B, Warmuth MK, editors. Learning Theory and Kernel Machines.
#' Berlin, Heidelberg: Springer Berlin Heidelberg; 2003. p. 144–58.
#' (Goos G, Hartmanis J, van Leeuwen J, editors. Lecture Notes in Computer
#' Science; vol. 2777). http://link.springer.com/10.1007/978-3-540-45167-9_12
#'
#' @return
#'
#' a list object containing the formula call and results by pathway
#'
#' @examples
#' data(smokers)
#'
#' # Query KEGGREST API
#' networks <- GatherNetworks(SE = smokers, keggID = "kegg_id",
#' species = "hsa", minPathwaySize = 5)
#'
#' # Run PaIRKAT Analysis
#' output <- PaIRKAT(log_FEV1_FVC_ratio ~ age, networks = networks)
#'
#' # View Results
#' output$results
#'
#' @export
#'
PaIRKAT <- function(formula.H0, networks, tau = 1) {

    # check if a valid networks object has been passed
    if(!is(networks, "list")){
        stop("Invalid `networks` object. Pass a valid object created with
             GatherNetworks function")
    }
    if(!all(c("networks","pdat","SE") %in% names(networks))){
        stop("Invalid `networks` object. Pass a valid object created with
             GatherNetworks function")
    }
    if(length(networks$networks) == 0){
        stop("`networks` object contains 0 pathways")
    }

    # Unpack SE object into phenotype, metabolite, and pathway data
    SE <- networks$SE
    mD <- assays(SE)[[1]]
    tmD <- transpose_tibble(mD)
    # pD <- tibble::as_tibble(rowData(SE), .name_repair = "minimal")
    # pD$rowname <- rownames(rowData(SE))
    # pD <- tibble::column_to_rownames(pD, var = "rowname")
    cD <- tibble::as_tibble(colData(SE), .name_repair = "minimal")
    cD$rowname <- rownames(colData(SE))
    cD <- tibble::column_to_rownames(cD, var = "rowname")

    # check formula.H0 for proper formatting
    if (!is(formula.H0, "formula")){
        stop("Invalid `formula.H0`. Please format
             as `outcome ~ covariate_1 + covariate_2 + ... covariate_n`")
    }

    # parse variables from formula
    keepVars <- paste(formula.H0[[2]])
    if (length(formula.H0[[3]]) > 1) {
        for (i in seq_len(length(formula.H0[[3]]))) {
            keepVars <- c(keepVars, paste(formula.H0[[3]][[i]]))
        }
    }
    else {
        keepVars <- c(keepVars, paste(formula.H0[[3]]))
    }

    # remove + from list of variables
    keepVars <- keepVars[keepVars != "+"]

    # check if variables are in the data
    if(!all(keepVars %in% names(cD))){
        stop(keepVars[!(keepVars %in% names(cD))]," not found in data")
    }

    cD <- cD[, names(cD) %in% keepVars]

    # check for missing data and subset complete cases
    completeCases <- complete.cases(cD)
    cD <- cD[completeCases,]
    tmD <- tmD[completeCases,]

    # check properties of outcome and format or throw appropriate errors
    outcome_vector <- cD[[keepVars[1]]]

    if (length(unique(outcome_vector)) == 2) {
        out.type <- "binary"
        message("Binary outcome detected.
                Null model will be a glm with logit link")
    } else if (length(unique(outcome_vector)) > 2) {
        out.type <- "continuous"
        message("Continuous outcome detected.
                Null model will be a linear model.")
    }

    if (!is.numeric(outcome_vector)) {
        if (out.type == "binary") {
            ref_level <- unique(outcome_vector)[1]
            pos_level <- unique(outcome_vector)[2]
            outcome_vector[outcome_vector == ref_level] <- 0
            outcome_vector[outcome_vector == pos_level] <- 1
            cD[keepVars[1]] <- outcome_vector

            warning(
                "Binary outcome is non-numeric, encoding",
                pos_level,
                "as 1 and",
                ref_level,
                "as 0."
            )
        }
        else{
            stop("Invalid formula: outcome is non-numeric")
        }

    }

    # Perform kernel tests
    pp_frame <- NULL

    for (i in seq_len(length(networks$networks))) {
        G <- networks$networks[[i]]
        varnames <- igraph::V(G)$label
        ZZ <- scale(tmD[, varnames[varnames %in% colnames(tmD)]])

        ## normalized Laplacian
        L <- igraph::graph.laplacian(G, normalized = TRUE)
        rho <- median(dist(ZZ))

        Z <- ZZ %*% solve(diag(nrow(L)) + tau * L)

        K <- Gaussian_kernel(rho, Z)

        if (out.type == "continuous") {
            pp <- SKAT.c(formula.H0, .data = cD, K = K)
        }

        if (out.type == "binary") {
            pp <- SKAT.b(formula.H0, .data = cD, K = K)
        }

        pp$pathway <- names(networks$networks[i])
        pp_frame <- rbind(pp_frame, as.data.frame(pp))


    }
    list(call = formula.H0, results = pp_frame[, c(3, 2, 1)])
}

# Helpers -----------------------------------------------------------------

## Function for making formula from "..." in functions
formula_fun <- function(Y, covs) {
    cc <- character(0)
    if (length(covs) > 1) {
        ff <- paste("~", paste(paste0("`", covs, "`"), collapse = "+"))
    } else if (length(covs) == 1) {
        ff <- paste0("~ `", covs, "`")
    } else{
        ff <- "~ 1"
    }
    paste(paste0("`", Y, "`"), ff)
}

Gaussian_kernel <- function(rho, Z) {
    exp(-(1 / rho) * as.matrix(dist(
        Z, method = "euclidean", upper = TRUE
    ) ^
        2))
}

# Davies Test -------------------------------------------------------------

# The following code is taken from:
# https://github.com/jchen1981/SSKAT/blob/main/R/SSKAT.R
# Manuscript can be found at
# https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.21934

#Compute the tail probability of 1-DF chi-square mixtures
#' @importFrom CompQuadForm davies
KAT.pval <- function(Q.all,
                     lambda,
                     acc = 1e-9,
                     lim = 1e6) {
    pval = rep(0, length(Q.all))
    i1 = which(is.finite(Q.all))
    for (i in i1) {
        tmp <- CompQuadForm::davies(Q.all[i], lambda, acc = acc, lim = lim)
        pval[i] = tmp$Qq

        if (tmp$ifault > 0)
            warning("ifault = ", tmp$ifault)
    }
    return(pval)
}

SKAT.c <- function(formula.H0,
                   .data = NULL,
                   K,
                   acc = 0.00001,
                   lim = 10000,
                   tol = 1e-10) {
    formula.H0 <- formula(formula.H0)
    m0 <- lm(formula.H0, data = .data)
    mX <- model.matrix(formula.H0, data = .data)

    res <- resid(m0)
    df <- nrow(mX) - ncol(mX)
    s2 <- sum(res ^ 2)

    P0  <- diag(nrow(mX)) - mX %*% (solve(t(mX) %*% mX) %*% t(mX))
    PKP <- P0 %*% K %*% P0
    q <- as.numeric(res %*% K %*% res / s2)

    ee <- eigen(PKP - q * P0, symmetric = TRUE)
    lambda <- ee$values[abs(ee$values) >= tol]

    p.value <- KAT.pval(
        0,
        lambda = sort(lambda, decreasing = TRUE),
        acc = acc,
        lim = lim
    )

    return(list(p.value = p.value, Q.adj = q))
}

SKAT.b <- function(formula.H0,
                   .data = NULL,
                   K,
                   acc = 0.00001,
                   lim = 10000,
                   tol = 1e-10) {
    formula.H0 <- formula(formula.H0)
    X1 <- model.matrix(formula.H0, .data)
    lhs <- formula.H0[[2]]
    y <- eval(lhs, .data)

    y <- factor(y)


    if (nlevels(y) != 2) {
        stop('The phenotype is not binary!\n')
    } else {
        y <- as.numeric(y) - 1
    }

    glmfit <- glm(y ~ X1 - 1, family = binomial)

    betas <- glmfit$coef
    mu  <- glmfit$fitted.values
    eta <- glmfit$linear.predictors
    res.wk <- glmfit$residuals
    res <- y - mu

    w   <- mu * (1 - mu)
    sqrtw <- sqrt(w)

    adj <- sum((sqrtw * res.wk) ^ 2)

    DX12 <- sqrtw * X1

    qrX <- qr(DX12)
    Q <- qr.Q(qrX)
    Q <- Q[, seq_len(qrX$rank), drop = FALSE]

    P0 <- diag(nrow(X1)) - Q %*% t(Q)

    DKD <- tcrossprod(sqrtw) * K
    tQK <- t(Q) %*% DKD
    QtQK <- Q %*% tQK
    PKP <- DKD - QtQK - t(QtQK) + Q %*% (tQK %*% Q) %*% t(Q)
    q <- as.numeric(res %*% K %*% res) / adj
    ee <- eigen(PKP - q * P0, symmetric = TRUE, only.values = TRUE)
    lambda <- ee$values[abs(ee$values) >= tol]

    p.value <- KAT.pval(
        0,
        lambda = sort(lambda, decreasing = TRUE),
        acc = acc,
        lim = lim
    )

    return(list(p.value = p.value, Q.adj = q))
}

# Saddle pVal functions ---------------------------------------------------

saddle = function(x, lambda) {
    d = max(lambda)
    lambda = lambda / d
    x = x / d
    k0 = function(zeta)
        - sum(log(1 - 2 * zeta * lambda)) / 2
    kprime0 = function(zeta)
        sapply(zeta,
               function(zz)
                   sum(lambda / (1 - 2 * zz * lambda)))
    kpprime0 = function(zeta)
        2 * sum(lambda ^ 2 / (1 - 2 * zeta * lambda) ^ 2)
    n = length(lambda)
    if (any(lambda < 0)) {
        lmin = max(1 / (2 * lambda[lambda < 0])) * 0.99999
    } else if (x > sum(lambda)) {
        lmin = -0.01
    } else {
        lmin = -length(lambda) / (2 * x)
    }
    lmax = min(1 / (2 * lambda[lambda > 0])) * 0.99999
    hatzeta = uniroot(
        function(zeta)
            kprime0(zeta) - x,
        lower = lmin,
        upper = lmax,
        tol = 1e-08
    )$root
    w = sign(hatzeta) * sqrt(2 * (hatzeta * x - k0(hatzeta)))
    v = hatzeta * sqrt(kpprime0(hatzeta))
    if (abs(hatzeta) < 1e-4) {
        return(NA)
    } else{
        return(pnorm(w + log(v / w) / w, lower.tail = FALSE))
    }
}

Liu.pval = function(Q, lambda) {
    c1 = rep(0, 4)
    for (i in seq_len(4)) {
        c1[i] = sum(lambda ^ i)
    }
    muQ = c1[1]
    sigmaQ = sqrt(2 * c1[2])
    s1 = c1[3] / c1[2] ^ (3 / 2)
    s2 = c1[4] / c1[2] ^ 2
    if (s1 ^ 2 > s2) {
        a = 1 / (s1 - sqrt(s1 ^ 2 - s2))
        d = s1 * a ^ 3 - a ^ 2
        l = a ^ 2 - 2 * d
    } else {
        l = 1 / s2
        a = sqrt(l)
        d = 0
    }
    muX = l + d
    sigmaX = sqrt(2) * a

    Q.Norm = (Q - muQ) / sigmaQ
    Q.Norm1 = Q.Norm * sigmaX + muX
    pchisq(Q.Norm1,
           df = l,
           ncp = d,
           lower.tail = FALSE)
}

Sadd.pval = function(Q.all, lambda) {
    sad = rep(1, length(Q.all))
    if (sum(Q.all > 0) > 0) {
        sad[Q.all > 0] = sapply(Q.all[Q.all > 0], saddle, lambda = lambda)
    }
    id = which(is.na(sad))
    if (length(id) > 0) {
        sad[id] = Liu.pval(Q.all[id], lambda)
    }
    return(sad)
}
