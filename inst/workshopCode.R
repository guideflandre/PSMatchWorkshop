### Setup and Installation ###

## Install latest versions from GitHub (for validatePSM and other new features)
remotes::install_github(c(
    "rformassspectrometry/PSMatch@validatePSM",
    "rformassspectrometry/PTMods"
), build_vignettes = TRUE, force = TRUE)

## Install other required Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(c("Spectra", "MsDataHub", "TargetDecoy"))


### Load Libraries ###

library(PSMatch)
library(PTMods)
library(Spectra)


### Load Example Data ###

## Load psmBoekweg, the identification file as a result from a database search
data("psmBoekweg")


### Create PSM Object ###

## Create a PSM/spectrum identifier
head(psmBoekweg$pkey <- paste0(
    basename(psmBoekweg$filename),
    sub("^.+scan=", "::", psmBoekweg$scannr)
))

psms <- PSM(
    x = psmBoekweg,
    spectrum = "pkey",
    peptide = "peptide",
    protein = "protein",
    decoy = "label",
    rank = "rank",
    score = "hyperscore",
    fdr = "spectrum_q"
)
psms

## Load the TMT Erwinia mzIdentML file
f <- MsDataHub::TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.20141210.mzid()
(
    id <- PSM(
        x = f,
        spectrum = "spectrumID",
        peptide = "sequence",
        protein = "DatabaseAccess",
        decoy = "isDecoy",
        rank = "rank",
        score = "MS.GF.RawScore",
        fdr = "MS.GF.QValue"
    )
)


### Explore PSM Data Structure ###

## Display PSM variable names used in the package
(psmVar <- psmVariables(psms))
(idVar <- psmVariables(id))
names(psms)


### Decoy and Target Hits ###

## Evaluate target-decoy distribution
library(TargetDecoy)
evalTargetDecoysHist(data.frame(psms),
    decoy = "label", score = "sage_discriminant_score",
    log10 = FALSE, nBins = 80
)

## Filter out decoy hits
psms
(psmsNoDecoys <- filterPsmDecoy(psms))


### Multiple Matches Per Spectrum ###

## Count how many matches per spectrum
table(table(id$spectrumID))

## Inspect specific scan with multiple matches
i <- grep("scan=1774", id$spectrumID)
data.frame(id)[i, c(idVar["spectrum"], "sequence", "DatabaseAccess", "rank")]


### Reduce PSMs to One Row Per Spectrum ###

## Collapse multiple matches into single row
dim(id)
idReduced <- reducePSMs(id, id$spectrumID)
dim(idReduced)

## Verify no duplicates
anyDuplicated(idReduced$spectrumID)

## Check list structure for multiple values
j <- grep("scan=1774", idReduced$spectrumID)
data.frame(idReduced)[j, "DatabaseAccess"]

## Check if PSM object is reduced
reduced(idReduced)


### Filter PSM data ###

## Apply filterPSMs and filterPsmFdr
id |> filterPSMs()
psms |> filterPsmFdr(FDR = 0.01)


### Describe peptides and proteins ###

describePeptides(psms)
describeProteins(psms)


### Adjacency matrices ###

## Reload and filter keeping shared peptides
psmsFiltered <- psms |>
    filterPsmDecoy() |>
    filterPsmRank()
data.frame(psmsFiltered[1:10, psmVar[c("peptide", "protein")]])

## Build adjacency matrix (weighted by score)
adj <- makeAdjacencyMatrix(psmsFiltered)
dim(adj)
adj[1:5, 1:5]

## Set score = NA for counts only and binary = TRUE for a binary matrix
adj <- makeAdjacencyMatrix(psmsFiltered, score = NA)
dim(adj)
adj[1:5, 1:5]


### Connected Components Analysis ###

## Decompose adjacency matrix into connected components
cc <- ConnectedComponents(adj)
cc
length(cc)

## Examples of different connected component types
connectedComponents(cc, 2)
connectedComponents(cc, 22)
connectedComponents(cc, 1)
connectedComponents(cc, 11)


### Visualise complex components ###

## Find large connected components
head(largeCC <- which(nrows(cc) > 2 & ncols(cc) > 2), 15)
dims(cc[head(largeCC, 15)])

## Extract and plot specific component
cx <- connectedComponents(cc, 1049)
cx
plotAdjacencyMatrix(cx)


### Prioritise connected components ###

## Compute metrics and rank components
cctab <- prioritiseConnectedComponents(cc)
head(cctab, 15)

## PCA visualisation of component metrics
library("factoextra")
fviz_pca(prcomp(cctab, scale = TRUE, center = TRUE))


### PTM annotation conversion ###

## Name → deltaMass
convertAnnotation(
    "M[Oxidation]PEPTIDE",
    convertToStyle = "deltaMass"
)

## deltaMass → name
convertAnnotation(
    "M[+15.995]PEPTIDE",
    convertToStyle = "name"
)


### Add fixed modifications ###

## Add carbamidomethylation to cysteines
addFixedModifications(
    "SCALITDCGR",
    fixedModifications = c(C = "Carbamidomethyl")
)

## TMT labelling of the N-terminus (mass 304.207 Da)
addFixedModifications(
    "SCALITDGR",
    fixedModifications = c(Nterm = 304.207, C = "Carbamidomethyl")
)

## Fluoro modification at position 1 only
addFixedModifications(
    "CALSITDGR",
    fixedModifications = c("Fluoro"),
    pos = 4
)


### Add variable modifications ###

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

data("spBoekweg")

head(spBoekweg$pkey <- paste0(
    basename(spBoekweg$dataOrigin),
    sub("^.+scan=", "::", spBoekweg$spectrumId)
))

## Merge by spectrum identifier
sp <- joinSpectraData(spBoekweg, psms,
    by.x = "pkey",
    by.y = "pkey"
)
sp

## Set 'sequence' variable equal to 'peptide' but with converted PTM notation
sp$sequence <- convertAnnotation(sp$peptide, "name")


### Select Single PSM for Analysis ###

## Filter to specific scan
sp6570 <- filterPrecursorScan(sp, 6570)
(sek <- sp6570$sequence)


### Calculate Theoretical Fragment Ions ###

## Compute fragment ions for peptide sequence
frags <- calculateFragments(
    sek[2],
    addCarbamidomethyl = FALSE
)
frags


### Calculate Fragments with PTM Enumeration ###

## Generate variable modification combinations
(
    varSeqs <- PTMods::addVariableModifications("ATK",
        variableModifications = c(T = "Phospho")
    )
)

calculateFragments(varSeqs, verbose = FALSE)


### Visualise Spectra ###

## Plot raw spectrum without annotation
plotSpectra(sp6570[2], allCharges = TRUE)

## Plot annotated spectrum with PTMs
plotSpectraPTM(sp6570[2], allCharges = TRUE)


### Adjust PPM Tolerance ###

## Increase ppm tolerance to match more peaks
plotSpectraPTM(sp6570[2], ppm = 150, allCharges = TRUE)


### Compare Modification Scenarios ###

## Compare spectra with and without oxidation
spComp <- sp6570[2]
spComp$sequence <- getCanonicalSequence(spComp$sequence)
plotSpectraPTM(spComp,
    variableModifications = c(M = "Oxidation"),
    deltaMz = FALSE,
    asp = 2 / 1,
    main = c(
        "Scan 6570, no oxidation",
        "Scan 6570, with oxidation"
    )
)


### Validate PSMs ###

scans_seq <- precScanNum(sp[!is.na(sp[["sequence"]])])
sp <- filterPrecursorScan(sp, scans_seq)

(
    res <- validatePSM(sp[1:20],
        peptideVariable = psmVar["peptide"],
        fdr = psmVar["fdr"]
    )
)

## Add metrics to initial PSM object
psm_df <- as.data.frame(psms)
combined <- merge(psm_df, res,
    by.x = "scannr", by.y = "spectrumId",
    all.x = FALSE
)

## Check overlap of matched fragments
sp_ms2 <- sp[msLevel(sp) == 2L & !is.na(sp$sequence)]

checkOverlap(sp_ms2[10])
plotSpectraPTM(sp_ms2[10])


### Session information ###

sessionInfo()
