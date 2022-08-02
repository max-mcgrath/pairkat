#' @title Gather pathway information from RaMP-DB
#' @name GatherNetworks2
#' @description Takes a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'  and constructs a list of KEGG or RaMP pathways along with their associated 
#'  networks
#'  
#' @param SE A \code{\link[SummarizedExperiment]{SummarizedExperiment}} with
#'     the features of interest in the 
#'     first \code{\link[SummarizedExperiment]{assay}}.
#' @param idCols Character vector indicating the name of the column(s) in \code{SE}
#'     row data containing database IDs
#' @param database Character vector indicating the name of the databases you 
#'     wish to query. Options are "KEGG", "RaMP", and "RaMPLocal"
#' @param minPathwaySize Restriction on the minimum number of nodes in networks
#'     returned. Must be greater than or equal to 2.
#' @param species The three letter KEGG organism ID. If using RaMP, this
#'     argument is ignored.
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
#' # Gather networks from RaMP-DB
#' networks <- GatherNetworksRaMPLocal(SE = smokers, idCols = "kegg_id",
#'                                     database = "KEGG",
#'                                     minPathwaySize = 5)
#' }
#' 
#' @importFrom KEGGREST keggList keggGet keggLink
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom tibble as_tibble column_to_rownames rownames_to_column
#' @importFrom stats binomial complete.cases dist formula glm lm median 
#'     model.matrix pchisq pnorm resid uniroot
#' @importFrom igraph V graph.laplacian graph_from_adjacency_matrix 
#'     graph_from_data_frame delete_vertices vertices gorder simplify
#' @importFrom data.table transpose
#' @importFrom methods is
#' @importFrom magrittr %>%
#' @importFrom dplyr select filter left_join distinct relocate all_equal
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @export
GatherNetworks2 <- function(SE,
                            idCols = "KEGG",
                            database = c("KEGG", "RaMP"),
                            nameCol = idCols[1],
                            species = "hsa",
                            minPathwaySize = 5) {
    # Argument preparation
    database <- match.arg(database)
    call <- as.list(match.call())
    call <- list(call = call[[1]],
                 idCols = idCols,
                 database = database,
                 nameCol = nameCol,
                 species = species,
                 minPathwaySize = minPathwaySize)
    
    # Input validation
    if (!is(SE, "SummarizedExperiment")) {
        stop("SE object is not a SummarizedExperiment")
    }
    
    if (as.integer(minPathwaySize) < 2) {
        stop("Minimum pathway size should be 2 or more")
    }
    if (!all(idCols %in% colnames(rowData(SE)))) {
        stop("ID(s): ", idCols[!(idCols %in% colnames(rowData(SE)))],
             " not found in SummarizedExperiment object rowData")
    }
    if (!(nameCol %in% colnames(rowData(SE)))) {
        warning("nameCol: ", nameCol, 
                " not found in SummarizedExperiment object rowData, ",
                " defaulting to first idCol (", idCols[1])
        nameCol <- idCols[1]
    }
    
    # Remove metabolites with zero variance
    SE <- .removeConstMets(SE)
    
    # Extract database IDs from SE, store in appropriate form (vector for KEGG,
    #   data frame for RaMP)
    metabolites <- .prepareMetabolites(SE, idCols, database)
    
    # Get networks from appropriate database(s)
    if(database == "KEGG") {
        keggRtn <- .getNetworksKEGG(metabolites, species, nameCol, 
                                    minPathwaySize)
        rampRtn <- NULL
    } else if (database == "RaMP") {
        rampRtn <- .getNetworksRaMP(metabolites, idCols, nameCol, 
                                    minPathwaySize)
        keggRtn <- NULL
    }
    
    # Return list networks, pdat, and SE
    list(networks = c(keggRtn$networks, rampRtn$networks), 
         pdat = c(keggRtn$pathData, rampRtn$pathData), SE = SE,
         call = call)
    
}

# Remove metabolites with zero variance
.removeConstMets <- function(SE) {
    
    # Determine indices to remove
    constIndices <- 
        which(apply(assay(SE), 1, function(x) length(unique(x)) == 1))
    
    # Throw warning, remove indices which are constant
    if (length(constIndices) > 0) {
        warning( "Removed ", length(constIndices), 
                 " compounds due to insufficient variance")
        SE <- SE[-constIndices, ]
    }
    
    SE
}


.prepareMetabolites <- function(SE, idCols, database) {
    
    metabolites <- as.data.frame(rowData(SE))
    
    # Reduce to metabolites which have at least one database ID
    indexNARows <- dplyr::select(metabolites, all_of(idCols)) |> 
        apply(1, function(row) !all(is.na(row)))
    metabolites <- metabolites[indexNARows, ]
    
    # Add internal uniqueID (essentially row number, but won't change as
    #   DF gets chopped up and rearranged)
    metabolites$internalID <- 1:nrow(metabolites)
    metabolites <- metabolites |> relocate(internalID)
}

.getNetworksKEGG <- function(metabolites, species, nameCol, minPathwaySize) {
    
    # Create data frame of reaction/compound pairs - used by .createKEGGNetworks
    #   to identify edges
    reactions <- keggLink("compound", "reaction")
    reactions <- substr(reactions, 5, nchar(reactions))
    reactions <- data.frame(reactionID = names(reactions),
                            metID = reactions)
    
    # Get all pathways
    pathwayIDs <- unique(keggLink("pathway", species))
    
    # Get pathway list from KEGG
    pathData <- .getKEGGPathways(pathwayIDs, 10)
    
    # Create networks from pathway data
    networks <- lapply(pathData, .createNetworks, metabolites, reactions,
                       minPathwaySize, nameCol)
    
    names(networks) <- unlist(lapply(pathData, function(x) x$NAME))
    
    # Remove NULL networks
    networks <- networks[!unlist(lapply(networks, is.null))]
    
    # Reduce PathData to those which remain in networks
    pathData <- pathData[unlist(lapply(pathData, function (x) { x$NAME })) %in% 
                             names(networks)]
    
    # Return list of networks, path data
    list(networks = networks,
         pathData = pathData)
    
}

.getNetworksRaMP <- function(metabolites, idCols, nameCol, minPathwaySize) {
    
    # Create data frame of reaction/compound pairs
    reactions <- .getReactionsRaMP(metabolites, idCols, 50)
    
    # Get list of pathways
    pathData <- .getRaMPPathways(metabolites, idCols, 50)
    
    # Create networks from pathway data
    networks <- lapply(pathData, .createNetworks, metabolites, reactions,
                       minPathwaySize, nameCol)
    
    # Label networks with name of pathway
    names(networks) <- names(pathData)
    
    # Remove NULL networks
    networks <- networks[!unlist(lapply(networks, is.null))]
    
    # Reduce PathData to those which remain in networks
    pathData <- pathData[names(pathData) %in% names(networks)]
    
    # Return list of networks, path data
    list(networks = networks,
         pathData = pathData)
    
}

