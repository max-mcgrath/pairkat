#' @title Gather pathway information from KEGG through the KEGGREST API.
#' @name GatherNetworks
#' @description Takes a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'  and constructs a list with KEGG pathways
#'
#' @param SE A \code{\link[SummarizedExperiment]{SummarizedExperiment}} with
#' the features of interest in the
#' first \code{\link[SummarizedExperiment]{assay}}.
#'
#' @param keggID column name in pathway data containing KEGG IDs
#' @param species The three letter KEGG organism ID
#' @param minPathwaySize Filter pathways that are below a minimum size
#'
#' @details
#' Queries KEGG database for known molecular interactions between
#' included metabolites via the KEGGREST API.
#'
#' @return
#' a list object containing the original SummarizedExperiment and
#' igraph network objects
#'
#' @examples
#'
#' library(SummarizedExperiment)
#' data(smokers)
#' # Query KEGGREST API
#'
#' networks <- GatherNetworks(SE = smokers, keggID = "kegg_id",
#' species = "hsa", minPathwaySize = 5)
#'
#' @export
#'
GatherNetworks <-
    function(SE,
             keggID = "KEGG",
             species = "hsa",
             minPathwaySize = 5) {
        # Input validation
        if (!is(SE, "SummarizedExperiment")) {
            stop("SE object is not a SummarizedExperiment")
        }
        
        valid_species <-
            as.data.frame(keggList("organism"))$organism
        if (!(species %in% valid_species)) {
            stop("Invalid species ID.
         Run `keggList('organism')` to view available options.")
        }
        
        if (as.integer(minPathwaySize) < 2) {
            stop("Minimum pathway size should be 2 or more")
        }
        
        # Unpack SE object into phenotype, metabolite, and pathway data
        mD <- assays(SE)[[1]]
        tmD <- transpose_tibble(mD)
        pD <-
            tibble::as_tibble(rowData(SE), .name_repair = "minimal")
        pD$rowname <- rownames(rowData(SE))
        pD <- column_to_rownames(pD, var = "rowname")
        cD <-
            tibble::as_tibble(colData(SE), .name_repair = "minimal")
        cD$rowname <- rownames(colData(SE))
        cD <- column_to_rownames(cD, var = "rowname")
        
        # Validate KEGG IDs
        if (!(keggID %in% colnames(pD))) {
            stop(keggID,
                 "not found in SummarizedExperiment object rowData")
        }
        
        # remove metabolites with zero variance
        remove_vars <-
            colnames(tmD)[!as.vector(unlist(lapply(tmD, function(x)
                length(unique(
                    x
                )) > 1)))]
        
        # throw warning if metabolites were removed
        if (length(remove_vars) > 0) {
            warning(
                "Removed the following metabolites from analysis due to
            insufficient variance:",
                paste(remove_vars, collapse = ", ")
            )
            tmD <- tmD[, !(colnames(tmD) %in% remove_vars)]
        }
        
        # compile pathway data
        pdat <- pathList(
            pathDat = pD,
            .pathID = keggID,
            ms = minPathwaySize,
            s_id = species
        )
        
        # gather networks
        nn <- getNetworks(
            pathDat = pD,
            metab = tmD,
            database = "KEGG",
            pdat = pdat,
            pathID = keggID
        )
        
        # add SE object to output
        nn$SE <- SE
        nn
    }

# Helpers -----------------------------------------------------------------

getNetworks <- function(pathDat,
                        metab,
                        database = "KEGG",
                        pdat,
                        pathID) {
    networks <- lapply(
        pdat$testPaths$keggPath,
        getNetwork,
        .comps = pdat$comps,
        .metab = metab,
        .pathDat = pathDat,
        compoundReaction = pdat$compoundReaction,
        .pathID = pathID
    )
    
    ## Naming list
    pathway_names <- !is.na(pdat$testPaths$pathwayNames)
    names(networks) <- pdat$testPaths$pathwayNames[pathway_names]
    
    ## Removing networks without connections (getNetwork returns -1)
    keepNet <- !sapply(networks, is.double)
    networks <- networks[keepNet]
    pdat$testPaths <- pdat$testPaths[keepNet,]
    pdat$comps <- pdat$comps[keepNet]
    
    return(list(networks = networks, pdat = pdat))
}

pathList <- function(pathDat, .pathID, ms, s_id) {
    .cr <- keggLink("compound", "reaction")
    hsapath <- unique(KEGGREST::keggLink("pathway", s_id))
    cr <- substr(.cr, 5, nchar(.cr))
    reactions <- names(cr)
    
    compId <- pathDat[, .pathID]
    compId <- unlist(strsplit(compId[!is.na(compId)], "[,]"))
    
    
    comps <-
        lapply(hsapath, function(p)
            try(keggGet(p)[[1]], silent = TRUE)
        )
    names(comps) <- hsapath
    keepPaths <- unlist(lapply(comps, function(p)
        is(p, "list")))
    comps <- comps[keepPaths]
    results <- data.frame(keggPath = names(comps),
                          stringsAsFactors  =  FALSE)
    comp <- lapply(comps, function(p)
        names(p$COMPOUND))
    compNames <- sapply(comps, function(p)
        p$NAME)
    results$inpathway <-
        sapply(comp, function(co)
            sum(compId %in% co))
    results$pathwayNames <- compNames
    
    testPaths <- results[results$inpathway >= ms,]
    
    ## making labels cleaner
    spec_look <- as.data.frame(keggList("organism")[, 2:3])
    sn <- spec_look[spec_look$organism == s_id, 2]
    testPaths$pathwayNames <- sub(paste(" -", sn), "",
                                  testPaths$pathwayNames,
                                  fixed = TRUE)
    
    pdat <- list(
        testPaths = testPaths,
        comps = comps[names(comps) %in% testPaths$keggPath],
        pathDat = pathDat[!is.na(pathDat[, .pathID]) &
                              !duplicated(pathDat[, .pathID]),],
        compoundReaction = cr
    )
    
    pdat
}

## Calculates Laplacian of metabolite pathway ## pathVar,
getNetwork <- function(pathId,
                       .comps,
                       .metab,
                       .pathDat,
                       compoundReaction,
                       .pathID) {
    target_compound <- names(.comps[[pathId]]$COMPOUND)
    
    cvnames <- .pathDat[.pathDat[[.pathID]] %in% target_compound,]
    
    varnames <- rownames(cvnames)
    
    path.v <- intersect(names(.metab), varnames)
    
    cnames <- cvnames[path.v, .pathID, drop = TRUE]
    
    ncomp <- length(cnames)
    reactions <- names(compoundReaction)
    
    A <- matrix(-1, length(cnames), length(cnames))
    diag(A) <- 0
    for (i in seq_len(ncomp)) {
        r1 <- reactions[which(compoundReaction == cnames[i])]
        j <- 1
        while (j < i) {
            r2 <- reactions[which(compoundReaction == cnames[j])]
            common <- intersect(r1, r2)
            
            if (length(common) > 0) {
                A[i, j] <- A[j, i] <- 1
            } else
                A[i, j] <- A[j, i] <- 0
            
            j <- j + 1
        }
    }
    
    if (sum(A) == 0)
        return(-1)
    G <- igraph::graph_from_adjacency_matrix(A, mode = "undirected")
    igraph::V(G)$label <- path.v
    return(G)
}

# custom function to transpose while preserving names
transpose_tibble <- function(df) {
    t_df <- data.table::transpose(df)
    colnames(t_df) <- rownames(df)
    rownames(t_df) <- colnames(df)
    t_df <- t_df %>%
        tibble::rownames_to_column() %>%
        tibble::as_tibble()
    t_df <- column_to_rownames(t_df, var = "rowname")
    return(t_df)
}
