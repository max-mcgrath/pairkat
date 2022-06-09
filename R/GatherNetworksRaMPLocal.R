#' @title Gather pathway information from RaMP-DB
#' @name GatherNetworksRaMPLocal
#' @description Takes a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'  and constructs a list of RaMP pathways along with their associated networks
#'  
#' @param SE A \code{\link[SummarizedExperiment]{SummarizedExperiment}} with
#'     the features of interest in the 
#'     first \code{\link[SummarizedExperiment]{assay}}.
#' @param ID Character vector indicating the name of the column(s) in \code{SE}
#'     row data containing database IDs
#' @param minPathwaySize Filter pathways that are below a minimum size
#'
#' @details
#' Queries local RaMP database for known molecular interactions between
#' included metabolites
#'
#' @return
#' a list object containing the original SummarizedExperiment and
#' igraph network objects
#'
#' @examples
#' \dontrun{
#' 
#' data(smokers)
#' 
#' # Create RaMP environment
#' pkg.globals <- setConnectionToRaMP(dbname = "ramp2", username = "root", 
#'                                    conpass = "", host = "localhost")
#' 
#' # Gather networks from RaMP-DB
#' networks <- GatherNetworksRaMPLocal(SE = smokers, ID = "kegg_id",
#'                                minPathwaySize = 5)
#' }
#' 
#' 
#' @importFrom SummarizedExperiment assays colData rowData
#' @importFrom tibble as_tibble column_to_rownames rownames_to_column
#' @importFrom stats binomial complete.cases dist formula glm lm median 
#'     model.matrix pchisq pnorm resid uniroot
#' @importFrom igraph V graph.laplacian graph_from_adjacency_matrix 
#'     graph_from_data_frame delete_vertices vertices gorder
#' @importFrom data.table transpose
#' @importFrom methods is
#' @importFrom magrittr %>%
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @export
GatherNetworksRaMPLocal <- function(SE, IDs, minPathwaySize = 5) {
    
    # TODO: Add check for acceptable prefixes (hmdb, kegg, etc.)
    # TODO: Figure out if RaMP has different species and adjust accordingly
    
    # Input validation
    if (!is(SE, "SummarizedExperiment")) {
        stop("SE object is not a SummarizedExperiment")
    }
    
    if (as.integer(minPathwaySize) < 2) {
        stop("Minimum pathway size should be 2 or more")
    }
    
    # Unpack SE object into phenotype, metabolite, and pathway data
    mD <- assays(SE)[[1]]
    tmD <- .transposeTibble(mD)
    pD <-
        tibble::as_tibble(rowData(SE), .name_repair = "minimal")
    pD$rowname <- rownames(rowData(SE))
    pD <- tibble::column_to_rownames(pD, var = "rowname")
    cD <-
        tibble::as_tibble(colData(SE), .name_repair = "minimal")
    cD$rowname <- rownames(colData(SE))
    cD <- tibble::column_to_rownames(cD, var = "rowname")
    
    # Validate Database IDs
    if (!all(IDs %in% colnames(pD))) {
        stop("ID(s): ",
             IDs[!(IDs %in% colnames(pD))],
             " not found in SummarizedExperiment object rowData")
    }
    
    # Remove metabolites with zero variance
    remove_vars <- colnames(tmD)[!as.vector(unlist(lapply(tmD, function(x)
        length(unique( x )) > 1)))]
    
    # Throw warning if metabolites were removed
    if (length(remove_vars) > 0) {
        warning(
            "Removed the following metabolites from analysis due to
            insufficient variance:",
            paste(remove_vars, collapse = ", ")
        )
        tmD <- tmD[, !(colnames(tmD) %in% remove_vars)]
    }
    
    # Create single vector of relevant IDs
    metabolites <- select(pD, all_of(IDs))
    metabolites <- unlist(metabolites)
    metabolites <- unname(metabolites)
    metabolites <- metabolites[!is.na(metabolites)]
    
    # Compile pathway data
    rtn <- .getNetworksLocal(metabolites = metabolites, tmD = tmD,
                        minPathwaySize = minPathwaySize)
    
    rtn$SE <- SE
    
    rtn
    
}

.getNetworksLocal <- function(metabolites, tmD, minPathwaySize) {
    
    # Create ramp connection (environment with parent namespace:RaMP)
    # pkg.globals <- setConnectionToRaMP(dbname = dbName, conpass = conPass)
    
    # Get pathways associated with analytes
    pathwaydfids <- getPathwayFromAnalyte(metabolites)
    
    # Create list of all pathway names
    pathways <- unique(pathwaydfids$pathwayName)
    
    # Get genes which catalyze analytes
    new.metabolites <- rampFastCata(analytes = metabolites)
    
    # Get edges - defined as analytes which are both catalyzed by the same gene
    #   and are within the same pathway
    a <- new.metabolites[, 1:2]
    b <- a
    names(b) <- c("Output_Analyte", "Input_CatalyzedBy_CommonName")
    c <- merge(a, b, all = T)
    d <- tolower(.colwiseLocal(c[c$Input_Analyte != c$Output_Analyte, 2:3]))
    filtered_pathways <- pathwaydfids[pathwaydfids$commonName %in% d,1:2]
    colnames(d) <- c("Input_Analyte", "Output_Analyte")
    # this part is slow, could be improved
    e <- merge(d, filtered_pathways, by.x = "Input_Analyte", 
               by.y = "commonName")
    f <- merge(e, filtered_pathways, by.x = "Output_Analyte", 
               by.y = "commonName")
    g <- f[f$pathwayName.x == f$pathwayName.y,] %>% na.omit()
    pathway_names <- unique(g$pathwayName.x)
    
    # Create network graphs for all pathways with at least 1 edge
    networks <- list()
    
    for (pathway in pathway_names){
        all_metabs_pathway <- pathwaydfids[pathwaydfids$pathwayName == pathway,]
        path_network_df <- g[g$pathwayName.x == pathway,]
        edges <- path_network_df[,1:2]
        
        edges <- .colwiseLocal(edges)
        path_graph <- graph_from_data_frame(edges, directed = FALSE)
        single_nodes <- 
            all_metabs_pathway$commonName[!(all_metabs_pathway$commonName %in% 
                                                path_network_df$Input_Analyte |
                                                all_metabs_pathway$commonName %in%  # This and next two lines added, need to check correctness
                                                path_network_df$Output_Analyte)]
        single_nodes <- unique(single_nodes)
        path_graph <- path_graph + vertices(single_nodes)
        
        # Remove edges vertices which do not appear in original data
        rmVertices <- names(V(path_graph))[!(names(V(path_graph)) %in% 
                                                 colnames(tmD))]
        path_graph <- delete_vertices(path_graph, rmVertices)
        
        
        networks[[pathway]] <- path_graph
    }
    
    # Remove pathways with fewer vertices than minPathwaySize
    networks <- lapply(networks, function (graph) {
        if (gorder(graph) < minPathwaySize) {
            rtn <- NULL
        } else {
            rtn <- graph
        }
        
        rtn
    })
    
    networks <- networks[lengths(networks) == 10]
    
    return(list(networks = networks))
}

# Function for getting unique pairs of metabolites from two columns
.colwiseLocal <- function(dat) { 
    unique(cbind(pmin(dat[,1], dat[,2]), pmax(dat[,1], dat[,2])))
}


# Custom function to transpose while preserving names
.transposeTibbleLocal <- function(df) {
    t_df <- data.table::transpose(df)
    colnames(t_df) <- rownames(df)
    rownames(t_df) <- colnames(df)
    t_df <- t_df %>%
        tibble::rownames_to_column() %>%
        tibble::as_tibble()
    t_df <- column_to_rownames(t_df, var = "rowname")
    return(t_df)
}