# Set up SE w/ first column no variance
SE <- SummarizedExperiment(assay = data.frame(s1 = c(1, 2, 5, 8),
                                              s2 = c(1, 3, 6, 8),
                                              s3 = c(1, 4, 7, 8)),
                           colData = data.frame(pheno1 = c("a", "b", "c"),
                                                pheno2 = c("d", "e", "f")),
                           rowData = data.frame(ID1 = c("id1", "id2", "id3",
                                                        "id4"),
                                                ID2 = c("id5", "id6", "id7",
                                                        "id8")))

# Unpack assay data
mD <- assays(SE)[[1]]
tmD <- pairkat:::.transposeTibble(mD)

# Test that assay was transposed correctly
test_that("Check .transposeTibble", {
    expect_true(all(rownames(tmD) == c("s1", "s2", "s3")))
    expect_true(all(is.null(colnames(tmD))))
    expect_true(all(tmD == matrix(c(1, 1, 1, 2, 3, 4, 5, 6, 7, 8, 8, 8),
                                  nrow = 4)))
})

# Test that first and last metabolite are removed (and no others)
test_that("Check .removeConstMets()", {
    expect_warning(SE <- .removeConstMets(SE))
    newSE <- SummarizedExperiment(assay = data.frame(s1 = c(2, 5),
                                                     s2 = c(3, 6),
                                                     s3 = c(4, 7)),
                                  colData = data.frame(pheno1 = c("a", "b", "c"),
                                                       pheno2 = c("d", "e", "f")),
                                  rowData = data.frame(ID1 = c("id2", "id3"),
                                                       ID2 = c("id6", "id7")))
    rownames(newSE) <- c(2, 3)
    expect_true(identical(rowData(SE), rowData(newSE)))
    expect_true(identical(colData(SE), colData(newSE)))
    expect_true(all(assay(SE) == assay(newSE)))
})