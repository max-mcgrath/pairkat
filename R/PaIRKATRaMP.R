#' @title Perform PaIRKAT on the output from the GatherNetworksRaMP function
#' @name PaIRKATRaMP
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
#' \dontrun{
#' data(smokers)
#' 
#' # Create RaMP environment
#' pkg.globals <- setConnectionToRaMP(dbname = "ramp2", username = "root", 
#'                                    conpass = "", host = "localhost")
#' 
#' # Gather networks from RaMP-DB
#' networks <- GatherNetworksRaMP(SE = smokers, ID = "kegg_id",
#'                                minPathwaySize = 5)
#'                                
#' # Run PaIRKAT Analysis
#' output <- PaIRKAT(log_FEV1_FVC_ratio ~ age, networks = networks)
#'
#' # View Results
#' output$results
#' }
#'
#' @export
#'
PaIRKATRaMP <- function(formula.H0, networks, tau = 1) {
    
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
    tmD <- .transposeTibble(mD)
    pD <- tibble::as_tibble(rowData(SE), .name_repair = "minimal")
    pD$rowname <- rownames(rowData(SE))
    pD <- tibble::column_to_rownames(pD, var = "rowname")
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
        
        # Determine which DB IDs are present in network
        idsInNetwork <- igraph::vertex_attr(G) |>
            as.data.frame() |> 
            select(all_of(networks$call$idCols)) |>
            unlist(use.names = FALSE)
        
        # Determine the indices of rowData associated with above DB IDs
        relevantIndices <- as.data.frame(rowData(SE)) |> 
            select(all_of(networks$call$idCols)) |> 
            apply(1, function(row) { 
            any(idsInNetwork %in% row) })
        
        ZZ <- scale(tmD[, relevantIndices])
        
        ## normalized Laplacian
        L <- igraph::graph.laplacian(G, normalized = TRUE)
        rho <- median(dist(ZZ))
        
        tryCatch(Z <- ZZ %*% solve(diag(nrow(L)) + tau * L),
                 error = function(condition) {
                     message(paste0("Error with graph ", i))
                     message("Original error")
                     message(condition)
                 })
        
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


    