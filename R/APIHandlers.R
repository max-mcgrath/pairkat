.getKEGGPathways <- function(pathwayIDs, batchSize) {
    
    # Empty list to store pathway data, partition paths into batches of 10
    pathways <- list()
    batches <- 
        split(pathwayIDs, as.integer((seq_along(pathwayIDs) - 1) / batchSize))
    
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