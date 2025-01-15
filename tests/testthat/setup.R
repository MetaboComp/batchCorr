set.seed(2024)

################# Drift correction #################
batchB <- getBatch(
    peakTable = PTfill, meta = meta,
    batch = meta$batch, select = "B"
)
batchF <- getBatch(
    peakTable = PTfill, meta = meta,
    batch = meta$batch, select = "F"
)
# Drift correction using QCs
BCorr <- correctDrift(
    peakTable = batchB$peakTable,
    injections = batchB$meta$inj,
    sampleGroups = batchB$meta$grp, QCID = "QC",
    G = seq(5, 35, by = 3), modelNames = c("VVE", "VEE"),
    report = FALSE
)
# More unbiased drift correction using QCs & external reference samples
FCorr <- correctDrift(
    peakTable = batchF$peakTable,
    injections = batchF$meta$inj,
    sampleGroups = batchF$meta$grp, QCID = "QC",
    RefID = "Ref", G = seq(5, 35, by = 3),
    modelNames = c("VVE", "VEE"),
    report = FALSE
)

############### Merge batches ################
mergedData <- mergeBatches(list(BCorr, FCorr))

####### correctDrift low-level walkthrough ########
meta_test <- data.frame(
    injections = batchB$meta$inj,
    sampleGroups = batchB$meta$grp
)
batchQC <- .getGroup(
    peakTable = batchB$peakTable, meta = meta_test,
    sampleGroup = batchB$meta$grp, select = "QC"
)
QCObject <- makeQCObject(
    peakTable = batchQC$peakTable,
    inj = batchQC$meta$inj
)
clustered <- clust(QCObject$inj, QCObject$Feats,
    modelNames = c("VVE", "VEE"),
    G = seq(5, 35, by = 10), report = FALSE
)
calculated_loess <- driftCalc(clustered,
    smoothFunc = "loess",
    spar = 0.2, report = FALSE
)
calculated_spline <- driftCalc(clustered,
    smoothFunc = "spline",
    spar = 0.2, report = FALSE
)
BatchObject <- makeBatchObject(
    peakTable = batchB$peakTable,
    inj = batchB$meta$inj, QCObject = QCObject
)
corrected_none <- driftCorr(
    QCDriftCalc = calculated_spline,
    CorrObj = BatchObject, report = FALSE
)
cleaned_none <- cleanVar(corrected_none, CVlimit = 0.3, report = FALSE)

########## alignBatches low-level walkthrough #########
peakIn <- peakInfo(PT = PTnofill, sep = "@", start = 3)

bFlag <- batchFlag(
    PTnofill = PTnofill, batch = meta$batch,
    sampleGroup = meta$grp, peakInfo = peakIn, NAhard = 0.8
)

aIQ <- alignIndex(batchflag = bFlag, grpType = "QC", mzdiff = 0.002, rtdiff = 15, report = FALSE)

bAlign <- batchAlign(batchflag = bFlag, alignindex = aIQ, peaktable_filled = PTfill, batch = meta$batch)

############# SummarizedExperiment object ##############
peaks <- SimpleList(t(PTnofill), t(PTfill))
sampleData <- meta
featureData <- peakInfo(PT = PTnofill, sep = "@", start = 3)
rownames(featureData) <- rownames(peaks[[1]])
se <- SummarizedExperiment(assays = peaks, colData = sampleData, 
                           rowData = featureData)
names(assays(se)) <- c("nofill", "fill")
