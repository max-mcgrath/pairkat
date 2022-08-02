#' @title Gather pathway information from RaMP-DB
#' @name GatherNetworks2
#' @description Takes a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'  and constructs a list of KEGG or RaMP pathways along with their associated 
#'  networks
#'  
#' @param SE A \code{\link[SummarizedExperiment]{SummarizedExperiment}} with
#'     the features of interest in the 
#'     first \code{\link[SummarizedExperiment]{assay}}.
#' @param IDs Character vector indicating the name of the column(s) in \code{SE}
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
#' networks <- GatherNetworksRaMPLocal(SE = smokers, IDs = "kegg_id",
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
#'     graph_from_data_frame delete_vertices vertices gorder
#' @importFrom data.table transpose
#' @importFrom methods is
#' @importFrom magrittr %>%
#' @importFrom dplyr select filter left_join distinct relocate all_equal
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @export
GatherNetworks2 <- function(SE,
                            IDs = "KEGG",
                            database = c("KEGG", "RaMP"),
                            species = "hsa",
                            minPathwaySize = 5) {

    database <- match.arg(database)
    
    # Input validation
    if (!is(SE, "SummarizedExperiment")) {
        stop("SE object is not a SummarizedExperiment")
    }
    
    if (as.integer(minPathwaySize) < 2) {
        stop("Minimum pathway size should be 2 or more")
    }
    if (!all(IDs %in% colnames(rowData(SE)))) {
        stop("ID(s): ", IDs[!(IDs %in% colnames(rowData(SE)))],
             " not found in SummarizedExperiment object rowData")
    }
    
    # Remove metabolites with zero variance
    SE <- .removeConstMets(SE)
    
    # Extract database IDs from SE, store in appropriate form (vector for KEGG,
    #   data frame for RaMP)
    metIDs <- .prepareMetabolites(SE, IDs, database)
    
    # Get networks from appropriate database(s)
    if(database == "KEGG") {
        keggRtn <- .getNetworksKEGG(metIDs, species, minPathwaySize)
        rampRtn <- NULL
    } else if (database == "RaMP") {
        rampRtn <- .getNetworksRaMP(metIDs, IDs, minPathwaySize)
        keggRtn <- NULL
    }
    
    # Return list networks, pdat, and SE
    list(networks = c(keggRtn$networks, rampRtn$networks), 
         pdat = c(keggRtn$pathData, rampRtn$pathData), SE = SE)
    
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

# Extract database IDs from SE, store in appropriate form (vector for KEGG,
#   data frame for RaMP)
.prepareMetabolites <- function(SE, IDs, database) {
    
    metIDs <- as.data.frame(rowData(SE))
    
    if(database == "KEGG") {
        # Create single vector of relevant IDs to use for network creation
        metIDs <- dplyr::select(metIDs, all_of(IDs))
        metIDs <- unlist(metIDs, use.names = FALSE)
        metIDs <- metIDs[!is.na(metIDs)]
    } else if (database == "RaMP") {
        # Reduce to metabolites which have at least one database ID
        indexNARows <- dplyr::select(metIDs, all_of(IDs)) |> 
            apply(1, function(row) !all(is.na(row)))
        metIDs <- metIDs[indexNARows, ]
        
        # Add internal uniqueID (essentially row number, but won't change as
        #   DF gets chopped up and rearranged)
        metIDs$internalID <- 1:nrow(metIDs)
        metIDs <- metIDs |> relocate(internalID)
        
    } else {
        stop(database, " not a supported database")
    }
}

.getNetworksKEGG <- function(metIDs, species, minPathwaySize) {
    
    # Create data frame of reaction/compound pairs - used by .createKEGGNetworks
    #   to identify edges
    reactions <- keggLink("compound", "reaction")
    reactions <- substr(reactions, 5, nchar(reactions))
    reactions <- data.frame(reactionID = names(reactions),
                            metID = reactions)
    
    # Get all pathways
    pathwayIDs <- unique(keggLink("pathway", species))
    
    # Get pathway list from KEGG
    pathData <- .getKEGGPathways(metIDs, pathwayIDs)
    
    # Create networks from pathway data
    networks <- lapply(pathData, .createKEGGNetworks, metIDs, reactions,
                       minPathwaySize)
    
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

.getKEGGPathways <- function(metIDs, pathwayIDs) {
    
    # Empty list to store pathway data, partition paths into batches of 10
    pathways <- list()
    batches <- split(pathwayIDs, as.integer((seq_along(pathwayIDs) - 1) / 10))
    
    # Get pathways in batches
    for (batch in batches) {
        newPaths <- keggGet(batch)
        names(newPaths) <- batch
        pathways <- append(pathways, newPaths)
    }
    
    # Remove failed pathways
    keepPaths <- unlist(lapply(pathways, function(p) is(p, "list")))
    pathways <- pathways[keepPaths]
    
}

# Generate network based on which compounds are in pathway and data set -
#   meant to be lapply'd to list of pathways so pathway is single pathway obj
.createKEGGNetworks <- function(pathway, metIDs, reactions, minPathwaySize) {
    
    # Get compound IDs that are present in pathway
    metIDsInPathway <- names(pathway$COMPOUND)
    
    # Reduce to only those that are in OG data set
    metIDsInPathway <- intersect(metIDsInPathway, metIDs)
    
    # Return NULL if fewer remaining nodes than minPathwaySize
    if (length(metIDsInPathway) < minPathwaySize) return(NULL)
    
    # Reduce reactions to relevant ones
    reactions <- reactions[reactions$metID %in% metIDsInPathway, ]
    
    # Create data frame of edges
    edges <- dplyr::left_join(reactions, reactions, by = "reactionID")
    edges <- dplyr::filter(edges, metID.x != metID.y)
    edges <- dplyr::select(edges, x = metID.x, y = metID.y)
    edges <- dplyr::distinct(edges)
    
    # Create iGraph network/graph object
    network <- igraph::graph_from_data_frame(edges, directed = FALSE, 
                                             vertices = metIDsInPathway)
    
    # Remove duplicate edges
    network <- igraph::simplify(network)
}

.createRaMPNetworks <- function(pathway, metIDs, reactions, minPathwaySize) {
    
    
    # Get compound IDs that are present in pathway
    metIDsInPathway <- names(pathway$COMPOUND)
    
    # Reduce to only those that are in OG data set
    metIDsInPathway <- intersect(metIDsInPathway, unlist(metIDs))
    metIDsIndex <- apply(metIDs, 1, function(row) { 
        any(metIDsInPathway %in% row) })
    metIDs <- metIDs[metIDsIndex, ]
    
    # Return NULL if fewer remaining nodes than minPathwaySize
    if (nrow(metIDs) < minPathwaySize) return(NULL)
    
    # Reduce reactions to relevant ones
    reactions <- reactions[reactions$metID %in% metIDsInPathway, ]
    
    # Create data frame of edges
    edges <- dplyr::left_join(reactions, reactions, by = "reactionID")
    
    # Create data frame using unique IDs - this step is necessary to ensure we 
    #   use all information (multiple IDs per metabolite), but also do not
    #   end up with duplicate network structures
    leftNodeID <- apply(edges, 1, function(row) {
        metIDs$internalID[which(row[["metID.x"]] == metIDs, arr.ind = TRUE)[1]]
    })
    rightNodeID <- apply(edges, 1, function(row) {
        metIDs$internalID[which(row[["metID.y"]] == metIDs, arr.ind = TRUE)[1]]
    })
    newEdges <- data.frame(leftNodeID, rightNodeID) %>%
        dplyr::distinct() %>%
        dplyr::filter(leftNodeID != rightNodeID)
    
    # Create iGraph network/graph object
    network <- igraph::graph_from_data_frame(newEdges, directed = FALSE, 
                                             vertices = metIDs)
    
    # Remove duplicate edges
    network <- igraph::simplify(network)
}

.getNetworksRaMP <- function(metIDs, IDs, minPathwaySize) {
    
    # Create data frame of reaction/compound pairs
    reactions <- .getReactionsRaMP(metIDs, IDs, 50)
    
    # Get list of pathways
    pathData <- .getRaMPPathways(metIDs, IDs, 50)
    
    # Create networks from pathway data
    networks <- lapply(pathData, .createRaMPNetworks, metIDs, reactions,
                       minPathwaySize)
    
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

.getReactionsRaMP <- function(metIDs, IDs, batchSize) {
    
    # Create vector of all database IDs
    metIDs <- dplyr::select(metIDs, all_of(IDs))
    metIDs <- unlist(metIDs)
    names(metIDs) <- NULL
    metIDs <- metIDs[!is.na(metIDs)]
    
    # Partition metabolites into batches
    batches <- split(metIDs, as.integer((seq_along(metIDs) - 1) / batchSize))
    
    # Get genes which catalyze analytes in batches
    reactions <- lapply(batches, function(batch) {
        metJSON <- jsonlite::toJSON(list(analyte = batch))
        response <- httr::POST(
            "https://ramp-api-alpha.ncats.io/api/common-reaction-analytes",
            body = metJSON)
        newReactions <- jsonlite::fromJSON(rawToChar(response$content))$data
        newReactions
    })
    names(reactions) <- NULL
    
    # Reduce to single data frame with reaction/metabolite pairs
    reactions <- dplyr::bind_rows(reactions)
    reactions <- dplyr::select(reactions, metID = input_analyte, 
                          reactionID = rxn_partner_common_name,
                          metCommonName = input_common_names)
    reactions <- dplyr::distinct(reactions)
    
}

.getRaMPPathways <- function(metIDs, IDs, batchSize) {
    
    # Create vector of all database IDs
    metIDs <- dplyr::select(metIDs, all_of(IDs))
    metIDs <- unlist(metIDs)
    names(metIDs) <- NULL
    metIDs <- metIDs[!is.na(metIDs)]
    
    # Partition metabolites into batches
    batches <- split(metIDs, as.integer((seq_along(metIDs) - 1) / batchSize))
    
    # Get pathways in batches
    pathways <- lapply(batches, function(batch){
        metJSON <- jsonlite::toJSON(list(analytes = batch))
        response <- httr::POST(
            "https://ramp-api-alpha.ncats.io/api/pathways-from-analytes",
            body = metJSON)
        pathways <- jsonlite::fromJSON(rawToChar(response$content))$data
    })
    names(pathways) <- NULL
    
    pathways <- dplyr::bind_rows(pathways)
    pathways <- dplyr::distinct(pathways)
    
    # Transform pathway DF into list similar to KEGG version
    pathData <- lapply(unique(pathways$pathwayName), function(NAME) {
        
        pathway <- dplyr::filter(pathways, pathwayName == NAME)
        rtn <- list("NAME" = unique(pathway$pathwayName),
                    "SOURCE" = unique(pathway$pathwaySource),
                    "ENTRY" = unique(pathway$pathwayId),
                    "COMPOUND" = pathway$commonName)
        
        names(rtn$COMPOUND) <- pathway$inputId
        
        rtn
        
    })
    
    names(pathData) <- unique(pathways$pathwayName)
    
    pathData
    
}