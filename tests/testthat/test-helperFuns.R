# Set up SE w/ first column no variance
SE <- SummarizedExperiment::SummarizedExperiment(
    assay = data.frame(s1 = c(1, 4, 7, 9, 6, 0),
                       s2 = c(2, 5, 8, 8, 5, 0),
                       s3 = c(3, 6, 9, 7, 4, 0)),
    colData = data.frame(pheno1 = c("a", "b", "c"),
                         pheno2 = c("d", "e", "f")),
    rowData = data.frame(id1 = 
                             c("id1-1", "id1-2", "id1-3", "id1-4", NA, "id1-0"),
                         id2 = 
                             c("id2-1", "id2-2", "id2-3", NA, NA, "id2-0")))

# Test that last metabolite is removed (and no others)
test_that("Check .removeConstMets()", {
    expect_warning(SE <- pairkat:::.removeConstMets(SE))
    newSE <- SummarizedExperiment::SummarizedExperiment(
        assay = data.frame(s1 = c(1, 4, 7, 9, 6),
                           s2 = c(2, 5, 8, 8, 5),
                           s3 = c(3, 6, 9, 7, 4)),
        colData = data.frame(pheno1 = c("a", "b", "c"),
                             pheno2 = c("d", "e", "f")),
        rowData = data.frame(id1 = 
                                 c("id1-1", "id1-2", "id1-3", "id1-4", NA),
                             id2 = 
                                 c("id2-1", "id2-2", "id2-3", NA, NA)))
    expect_true(identical(SummarizedExperiment::rowData(SE), 
                          SummarizedExperiment::rowData(newSE)))
    expect_true(identical(SummarizedExperiment::colData(SE), 
                          SummarizedExperiment::colData(newSE)))
    expect_true(all(SummarizedExperiment::assay(SE) == 
                        SummarizedExperiment::assay(newSE)))
})

suppressWarnings(SE <- pairkat:::.removeConstMets(SE))
metIDsKEGG <- pairkat:::.prepareMetabolites(SE, c("id2"))
metIDsRaMP <- pairkat:::.prepareMetabolites(SE, c("id1", "id2"))

# Test that metabolite database IDs are properly prepared
test_that("Check .prepareMetabolites()", {
    expect_true(dplyr::all_equal(metIDsKEGG,
                                 data.frame(
                                     internalID = 1:3,
                                     id1 = c("id1-1", "id1-2", "id1-3"),
                                     id2 = c("id2-1", "id2-2", "id2-3"))))
    expect_true(dplyr::all_equal(
        metIDsRaMP,
        data.frame(internalID = 1:4,
                   id1 = c("id1-1", "id1-2", "id1-3", "id1-4"),
                   id2 = c("id2-1", "id2-2", "id2-3", NA))))
})

# Create test data corresponding to a simple network m3-m1-m2. However, two of
#   the nodes have multiple IDs. Functions should be able to handle this 
#   situation
pathway <- list(NAME = "testPathway",
                SOURCE = c("source1", "source2"),
                ENTRY = c("entryCode"),
                COMPOUND = c("met1", "met1", "met2", "met2", "met3"))
names(pathway$COMPOUND) <- c("id1-1", "id1-2", "id2-1", "id2-2", "id1-3")
metIDs <- data.frame(internalID = 1:3, 
                     db1 = c("id1-1", "id1-2", "id1-3"), 
                     db2 = c("id2-1", "id2-2", NA))
reactions <- data.frame(
    metID = c("id1-1", "id1-2", "id2-1", "id2-2", "id1-3", "id1-1"),
    reactionID = c(rep("reaction1", 4), rep("reaction2", 2)))
minPathwaySize <- 0

testNetwork <- pairkat:::.createNetworks(pathway, metIDs, reactions,
                                         minPathwaySize, "db1")

# Ensure that graphs have correct number of nodes, edges, and labels 
#   (attributes)
test_that("Check .createRaMPNetworks()", {
    expect_true(all(igraph::get.vertex.attribute(testNetwork, "db1") == 
                        c("id1-1", "id1-2", "id1-3")))
    expect_true(all(igraph::get.vertex.attribute(testNetwork, "db2")[1:2] == 
                        c("id2-1", "id2-2")))
    expect_true(all(igraph::gsize(testNetwork) == 2))
    expect_true(all(igraph::gorder(testNetwork) == 3))
})
