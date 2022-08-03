# Generate network based on which compounds are in pathway and data set -
#   meant to be lapply'd to list of pathways i.e. pathway is single pathway obj
.createNetworks <- function(pathway, metabolites, reactions, minPathwaySize,
                            nameCol) {
    
    # Get compound IDs that are present in pathway
    metabolitesInPathway <- names(pathway$COMPOUND)
    
    # Reduce to only those that are in OG data set
    metabolitesInPathway <- intersect(metabolitesInPathway, unlist(metabolites))
    metabolitesIndex <- apply(metabolites, 1, function(row) { 
        any(metabolitesInPathway %in% row) })
    metabolites <- metabolites[metabolitesIndex, ]
    
    # Return NULL if fewer remaining nodes than minPathwaySize
    if (nrow(metabolites) < minPathwaySize) return(NULL)
    
    # Reduce reactions to relevant ones
    reactions <- reactions[reactions$metID %in% metabolitesInPathway, ]
    
    # Create data frame of edges
    edges <- dplyr::left_join(reactions, reactions, by = "reactionID")
    
    # Create data frame using unique IDs - this step is necessary to ensure we 
    #   use all information (multiple IDs per metabolite), but also do not
    #   end up with duplicate network structures
    leftNodeID <- apply(edges, 1, function(row) {
        metabolites$internalID[which(row[["metID.x"]] == 
                                         metabolites, arr.ind = TRUE)[1]]
    })
    rightNodeID <- apply(edges, 1, function(row) {
        metabolites$internalID[which(row[["metID.y"]] == 
                                         metabolites, arr.ind = TRUE)[1]]
    })
    newEdges <- data.frame(leftNodeID, rightNodeID) %>%
        dplyr::distinct() %>%
        dplyr::filter(leftNodeID != rightNodeID)
    
    # Create iGraph network/graph object
    network <- igraph::graph_from_data_frame(newEdges, directed = FALSE, 
                                             vertices = metabolites)
    
    # Remove duplicate edges
    network <- igraph::simplify(network)
    
    # Update name attribute of igraph object
    names <- metabolites |> 
        select(all_of(nameCol)) |> 
        unlist(use.names = FALSE)
    igraph::V(network)$name <- names
    
    network
    
}