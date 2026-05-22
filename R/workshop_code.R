### Setup and Installation ###

## Install required Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("PSMatch", "PTMods", "Spectra", "MsDataHub", "TargetDecoy"))


### Load Libraries ###

library(PSMatch)
library(PTMods)
library(Spectra)
library(TargetDecoy)
library(ggplot2)


### Load Example Data ###

## Get example mzIdentML file from MsDataHub
f <- MsDataHub::TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.20141210.mzid()


### Create PSM Object ###

## Parse mzid file with PSM constructor
id <- PSM(f)

## Construct PSM object with specified column names
id <- PSM(x = f,
    spectrum = "spectrumID",
    peptide = "sequence",
    protein = "DatabaseAccess",
    decoy = "isDecoy",
    rank = "rank",
    score = "MS.GF.RawScore",
    fdr = "MS.GF.QValue")

id

### Explore PSM Data Structure ###

## Check dimensions of PSM object
dim(id)

## Display PSM variable names used in package
psmVar <- psmVariables(id)

## Show all column names
names(id)


### Decoy and Target Hits ###

## Count decoy vs target hits
table(id[[psmVar["decoy"]]])

## Evaluate target-decoy distribution
evalTargetDecoysHist(data.frame(id),
    decoy = "isDecoy", score = "MS.GF.RawScore",
    log10 = FALSE, nBins = 80
)


### Multiple Matches Per Spectrum ###

## Count how many matches per spectrum
table(table(id$spectrumID))

## Inspect specific scan with multiple matches
i <- grep("scan=1774", id$spectrumID)
data.frame(id)[i, c("sequence", "DatabaseAccess", "rank")]


### Reduce PSMs to One Row Per Spectrum ###

## Collapse multiple matches into single row
idReduced <- reducePSMs(id, id$spectrumID)
dim(idReduced)

## Verify no duplicates
anyDuplicated(idReduced$spectrumID)

## Check list structure for multiple values
j <- grep("scan=1774", idReduced$spectrumID)
data.frame(idReduced)[j, "DatabaseAccess"]

## Check if PSM object is reduced
reduced(idReduced)


### Filter PSM Data ###

## Remove decoy hits
id <- filterPsmDecoy(id)
nrow(id)
table(id$isDecoy)

## Keep only rank 1 matches
id <- filterPsmRank(id)
nrow(id)

## Example of shared peptide
id[id$sequence == "QKTRCATRAFKANKGRAR", "DatabaseAccess"]

## Remove shared peptides
id <- filterPsmShared(id)
nrow(id)

## Filter by FDR threshold
id <- filterPsmFdr(id)
nrow(id)


### Apply All Standard Filters at Once ###

## Reload and apply first three filters in one call
idFiltered <- PSM(f) |> filterPSMs()
nrow(idFiltered)


### Adjacency Matrices ###

## Reload and filter keeping shared peptides
id <- PSM(f) |> filterPsmDecoy() |> filterPsmRank()
data.frame(id[1:10, c("sequence", "DatabaseAccess")])

## Build adjacency matrix
adj <- makeAdjacencyMatrix(id)
dim(adj)
adj[1:5, 1:5]


### Connected Components Analysis ###

## Decompose adjacency matrix into connected components
cc <- ConnectedComponents(adj)
cc
length(cc)

## Examples of different connected component types
connectedComponents(cc, 1)
connectedComponents(cc, 38)
connectedComponents(cc, 527)
connectedComponents(cc, 920)


### Visualise Complex Components ###

## Find large connected components
largeCC <- which(nrows(cc) > 2 & ncols(cc) > 2)
dims(cc[largeCC])

## Extract and plot specific component
cx <- connectedComponents(cc, 1082)
cx
plotAdjacencyMatrix(cx)


### Prioritise Connected Components ###

## Compute metrics and rank components
cctab <- prioritiseConnectedComponents(cc)
head(cctab)

## PCA visualisation of component metrics
library("factoextra")
fviz_pca(prcomp(cctab, scale = TRUE, center = TRUE))


### PTM Annotation Conversion ###

## Convert between PTM annotation styles
convertAnnotation(
    "M[Oxidation]PEPTIDE",
    convertToStyle = "deltaMass"
)

convertAnnotation(
    "M[+15.995]PEPTIDE",
    convertToStyle = "name"
)


### Add Fixed Modifications ###

## Add carbamidomethylation to cysteines
addFixedModifications(
    "SCALITDGR",
    fixedModifications = c(C = "Carbamidomethyl")
)

## Add TMT labelling to N-terminus
addFixedModifications(
    "SCALITDGR",
    fixedModifications = c(Nterm = 304.207)
)

## Add modification at specific position
addFixedModifications(
    "SCALITDGR",
    fixedModifications = c("Fluoro"),
    pos = 1
)


### Add Variable Modifications ###

## Generate all possible modification combinations
addVariableModifications(
    "MPEPTIDE",
    variableModifications = c(M = -15.995, T = "Phospho")
)

## Limit maximum number of simultaneous modifications
addVariableModifications(
    "MPEPTIDE",
    variableModifications = c(M = -15.995, T = "Phospho"),
    maxMods = 1
)


### Strip Modifications ###

## Get canonical sequence without modifications
getCanonicalSequence("[+304]-M[Oxidation]PEPT[UNIMOD:21]IDE")


### Load Raw Spectra and Identification Data ###

## Load mzML file with raw spectra
spf <- MsDataHub::TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.20141210.mzML.gz()
sp  <- Spectra(spf)

## Load and filter identification data
idf <- MsDataHub::TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.20141210.mzid()
id  <- PSM(idf) |> filterPSMs()

## Merge spectra with identification data
sp <- joinSpectraData(sp, id,
    by.x = "spectrumId",
    by.y = "spectrumID"
)


### Select Single PSM for Analysis ###

## Filter to specific scan
sp1158 <- filterPrecursorScan(sp, 1158)
sp1158$sequence


### Calculate Theoretical Fragment Ions ###

## Compute fragment ions for peptide sequence
frags <- calculateFragments(
    sp1158$sequence[2],
    addCarbamidomethyl = FALSE
)
head(frags)


### Calculate Fragments with PTM Enumeration ###

## Generate variable modification combinations
varSeqs <- PTMods::addVariableModifications(
    sp1158$sequence[2],
    variableModifications = c(C = 57.02146, T = 79.966)
)
varSeqs

## Calculate fragments for all modified sequences
fragsVar <- calculateFragments(varSeqs)
nrow(fragsVar)
unique(fragsVar$peptide)


### Visualise Spectra ###

## Plot raw spectrum without annotation
plotSpectra(sp1158[2])

## Plot annotated spectrum with PTMs
dataOrigin(sp1158)[2] <- "TMT_Erwinia"
plotSpectraPTM(sp1158[2],
    fixedModifications = c(C = "Carbamidomethyl")
)


### Adjust PPM Tolerance ###

## Increase ppm tolerance to match more peaks
plotSpectraPTM(sp1158[2],
    fixedModifications = c(C = "Carbamidomethyl"),
    ppm = 50
)


### Compare Modification Scenarios ###

## Compare spectra with and without modification
plotSpectraPTM(sp1158[2],
    variableModifications = c(C = "Carbamidomethyl"),
    addCarbamidomethyl    = FALSE,
    deltaMz = FALSE,
    asp = 1/2,
    main = c(
        "Scan 1158 — no carbamidomethyl",
        "Scan 1158 — with carbamidomethyl"
    )
)


### Session Information ###

sessionInfo()
