test_that("Results from alignBatches has right dimensions", {
    alignBat <- alignBatches(
        peakInfo = peakIn, PeakTabNoFill = PTnofill,
        PeakTabFilled = PTfill, batches = meta$batch,
        sampleGroups = meta$grp, selectGroup = "QC",
        report = FALSE
    )

    # Check that all samples are preserved
    expect_identical(rownames(alignBat$PTalign), rownames(PTfill))
    # Check that there are less features after alignment
    expect_lt(ncol(alignBat$PTalign), ncol(PTfill))
    # Check that alignBatches output is of right length
    expect_equal(length(alignBat), 6)
})

test_that("alignBatches throws error when supplied various values", {
    # When there are missing values in PeakTabFilled
    expect_error(alignBatches(
        peakInfo = peakIn, PeakTabNoFill = PTnofill,
        PeakTabFilled = PTnofill, batches = meta$batch,
        sampleGroups = meta$grp, selectGroup = "QC",
        report = FALSE
    ))
    # When selectGroup value is not available in metadata
    expect_error(alignBatches(
        peakInfo = peakIn, PeakTabNoFill = PTnofill,
        PeakTabFilled = PTfill, batches = meta$batch,
        sampleGroups = meta$grp, selectGroup = "qc",
        mzdiff = -0.002, rtdiff = -2, report = FALSE
    ))
    # When sampleGroups == batches
    expect_error(alignBatches(
        peakInfo = peakIn, PeakTabNoFill = PTnofill,
        PeakTabFilled = PTfill, batches = meta$batch,
        sampleGroups = meta$batch, selectGroup = "QC",
        report = FALSE
    ))
    # When selectGroup value not available in metadata
    expect_error(alignBatches(
        peakInfo = peakIn, PeakTabNoFill = PTnofill,
        PeakTabFilled = PTfill, batches = meta$batch,
        sampleGroups = meta$grp, selectGroup = "sample"
    ))
})

test_that("alignBatches throws error when there are no alignment candidates", {
    # When PeakTabNoFill includes missing values
    expect_error(alignBatches(
        peakInfo = peakIn, PeakTabNoFill = PTfill,
        PeakTabFilled = PTfill, batches = meta$batch,
        sampleGroups = meta$grp, selectGroup = "QC",
        report = FALSE
    ))
    # When mzdiff and rtdiff are zero
    expect_error(alignBatches(
        peakInfo = peakIn, PeakTabNoFill = PTnofill,
        PeakTabFilled = PTfill, batches = meta$batch,
        sampleGroups = meta$grp, selectGroup = "QC",
        mzdiff = 0, rtdiff = 0, report = FALSE
    ))
    # When mzdiff and rtdiff are negative
    expect_error(alignBatches(
        peakInfo = peakIn, PeakTabNoFill = PTnofill,
        PeakTabFilled = PTfill, batches = meta$batch,
        sampleGroups = meta$grp, selectGroup = "QC",
        mzdiff = -0.002, rtdiff = -2, report = FALSE
    ))
    # When batches == sampleGroups
    expect_error(alignBatches(
        peakInfo = peakIn, PeakTabNoFill = PTnofill,
        PeakTabFilled = PTfill, batches = meta$grp,
        sampleGroups = meta$grp, selectGroup = "QC",
        report = FALSE
    ))
})

test_that("Output of batchFlag is correct", {
    # Check that object is right length
    expect_equal(length(bFlag), 4)

    bFlag_test <- cbind("batch" = c("B", "B", "F", "F", "H", "H"), "sampleGroup" = c("QC", "Ref", "QC", "Ref", "QC", "Ref"))
    # Check that batch and group data are correct
    expect_identical(bFlag$meta, bFlag_test)

    expected <- matrix(c(0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1),
        ncol = 3
    )
    # Check that flagged features are correct
    expect_identical(unname(bFlag$flagHard[1:6, 1:3]), expected)
})

test_that("Output of alignIndex is correct", {
    expected <- data.frame(
        "list" = c(0, 0, 3, 3, 0, 0),
        "sampleGroup" = c("", "", "QC", "QC", "", "")
    )
    # Check list of features to be combined
    expect_identical(unname(aIQ$shift[1:6, ]), unname(expected))

    expected <- data.frame(
        m = c(3, 73, 85, 119, 119, 120),
        n = c(4, 74, 87, 120, 121, 121)
    )
    # Check that misalignment candidates (m) and cluster candidates (n) are right
    expect_equal(unname(aIQ$events[1:6, c("m", "n")]), unname(expected))
    # Check the index of first cluster with three features
    expect_equal(which(aIQ$clusters$nFeat == 3)[1], 4)
})

test_that("Output of batchAlign is correct", {
    # Check index of first feature aligned with feature averaging
    expect_equal(which(bAlign$boolAveragedFill == TRUE)[1], 129)
    # Check index of first feature feature not kept after alignment
    expect_equal(which(bAlign$boolKeep == FALSE)[1], 4)
})

test_that("Assay control works (with alignBatches)", {
    # Assay not found
    expect_error(alignBatches(PeakTabNoFill = se, PeakTabFilled = se,
                     batches = "batch", sampleGroups = "grp", report = FALSE,
                     assay.type1 = "nope", assay.type2 = "fill", 
                     name = "aligned", rt_col = "rt", mz_col = "mz"))
    
    # Resultant assay must not be same as assay.type
    expect_error(alignBatches(PeakTabNoFill = se, PeakTabFilled = se,
                     batches = "batch", sampleGroups = "grp", report = FALSE,
                     assay.type1 = "nofill", assay.type2 = "fill", 
                     name = "fill", rt_col = "rt", mz_col = "mz"))
             
    # colData column not found
    expect_error(alignBatches(PeakTabNoFill = se, PeakTabFilled = se,
                     batches = "nope", sampleGroups = "grp", report = FALSE,
                     assay.type1 = "nofill", assay.type2 = "fill", 
                     name = "aligned", rt_col = "rt", mz_col = "mz"))
                 
    # rowData column not found
    expect_error(alignBatches(PeakTabNoFill = se, PeakTabFilled = se,
                     batches = "batch", sampleGroups = "grp", report = FALSE,
                     assay.type1 = "nofill", assay.type2 = "fill", 
                     name = "aligned", rt_col = "nope", mz_col = "mz"))
                     
    # Works with a single-assay objects without specifying assay.type nor name
    nofill <- se
    assays(nofill)[names(assays(nofill)) != "nofill"] <- NULL
    fill <- se
    assays(fill)[names(assays(fill)) != "fill"] <- NULL
    expect_no_error(alignBatches(PeakTabNoFill = nofill, PeakTabFilled = fill,
                        batches = "batch", sampleGroups = "grp", report = FALSE,
                        rt_col = "rt", mz_col = "mz"))
})
