# Handling proteomics identification data with PSMatch.

Authors: Guillaume Deflandre, Laurent Gatto and Sebastian Gibb.

## Overview 

### Description 

This workshop introduces you to the `PSMatch` _R_ and Bioconductor
package. As you will see, the package also depends on another Bioconductor
package called `PTMods` and we also highly suggest using the 
`Spectra` package for plotting spectra. All three packages are part
of the [R for Mass Spectrometry](https://www.rformassspectrometry.org/) package
series and you are welcome to learn more about this project by visiting its
corresponding website.

In short, `PSMatch` is a Bioconductor package for handling proteomics
identification data whilst `PTMods` handles post-translational modifications
(PTMs). You are encouraged to run the code interactively as I will present each
section.

By the end of this workshop, you will be able to import, explore, filter, and
visualise peptide-spectrum match (PSM) data; model peptide-protein
relationships using adjacency matrices; and handle PTM annotations in a
principled, reproducible manner.

### Pre-requisites 

- A conceptual understanding of shotgun proteomics: what a peptide
  is, what a mass spectrum is, and what a database search does. 

### Suggested readings 

- The `PSMatch` package vignettes, accessible via
`browseVignettes("PSMatch")`. 
- The pre-print available on the OSF platform:
[PSMatch: an R/Bioconductor package to explore proteomics identification
data.](https://doi.org/10.31219/osf.io/62v9p_v3). Note that the package has
gone through quite a few changes since.

### Participation 

As I already mentioned, you are welcome to follow along with me. Code is
provided in full, no typing from scratch is required. The code can be fetched
from [this GitHub repo](https://github.com/guideflandre/PSMatchWorkshop).

### Workshop goals 

- Understand the structure of peptide-spectrum match data and how it is
  represented in R using the `PSM` class.
- Recognise the problem of shared peptides and protein ambiguity, and
  know how adjacency matrices and connected components help address it.
- Appreciate the diversity of PTM annotation conventions and know how
  to convert between them and enumerate modified sequences.
- Understand the sometimes lacking confidence in identifications and as such,
  the importance of validating PSMs.

### Learning objectives 

- Import an identification data file into R and inspect the resulting `PSM`
  object.
- Apply the standard filtering steps (decoy removal, rank filtering,
  shared-peptide filtering, fdr filtering) individually and via `filterPSMs()`.
- Build an adjacency matrix and connected components from PSM data
  and visualise a protein group as a bipartite graph.
- Convert a PTM annotation between *deltaMass*, *unimodId*, and *name*
  formats using `convertAnnotation()`.
- Apply fixed and variable modifications to a peptide sequence and
  calculate theoretical fragment ions that account for those
  modifications.
- Plot spectra by including annotations based on the identification.
- Improve the confidence in a PSM using `validatePSM()`, built based
  on spectral and identification metrics.


## License 

The content of this workshop is provided under a
[CC-BY ShareAlike](https://creativecommons.org/licenses/by-sa/2.0/)
license.

To cite the PSMatch package in publications use:

```
  G. Deflandre, S. Gibbs and L. Gatto. PSMatch: an R/Bioconductor package to explore proteomics identification data.
 OSF;
  doi: https://doi.org/10.31219/osf.io/62v9p_v3

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {PSMatch: an R/Bioconductor package to explore proteomics identification data.},
    author = {Guillaume Deflandre and Sebastian Gibbs and Laurent Gatto},
    institution = {OSF},
    year = {2025},
    doi = {https://doi.org/10.31219/osf.io/62v9p_v2},
  }
```
