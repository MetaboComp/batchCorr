# batchCorr
Within and between batch correction of LC-MS metabolomics data

## Installation
Install `devtools` to be able to install packages from GitHub.

Install `batchCorr` package by:

`devtools::install_git("https://username:password@gitlab.com/CarlBrunius/batchCorr.git")`

where `username` and `password` are you actual username and password for your GitLab account.

In addition to functions relevant for within/between batch correction, data is provided to accurately reproduce figures from the original *Brunius et al* paper (see below).

## Workflow
After installation, a `Workflow_Example` folder is created in the `batchCorr` library (in your R library folder). Within this folder, there is a `workflow.R` script containing code on how the package was used to perform within/between batch correction and reproduce figures from the original *Brunius et al* paper (see below).

## Description
This is a repository containing functions within three areas of batch correction. These algorithms were originally developed 
to increase quality and information content in data from LC-MS metabolomics. However, the algorithms should be applicable to 
other data structures/origins, where within and between batch irregularities occur.

The three areas indicated are:

correction | abbreviation | description
:--- | :----------- | :----------
Batch alignment | BA | Functions to align features that are originally systematically misaligned between batches
Drift correction | DC | Functions to perform within batch intensity drift correction
Batch normalisation | BN | Funtions to perform between batch normalisation

### Batch alignment 
Batch alignment is achieved based on three concepts:
- Aggregation of feature presence/missingness on batch level.
- Identifying features with missingness within "the box", i.e. sufficiently similar in retention time and m/z.
- Ensuring orthogonal batch presence among feature alignment candidates.

### Drift correction
Drift correction is achieved based on:
- Clustering is performed on features in observation space (as opposed to the normally used observations in feature space)
- Clustering provides a tradeoff between 
  - modelling detail (multiple drift patterns within data set)
  - power per drift pattern
- Unbiased clustering is achieved using the Bayesian `mclust` R package

### Batch normalisation
Batch normalisation is achieved based on:
- QC/Reference (standard normalisation) or
- Population (median normalisation)
- The choice between the two is based on a quality heuristic determining whether the QC/Ref is suitable for normalisation. Otherwise population normalisation is performed instead.

## Reference
The development and inner workings of these algorithms are reported in:

*Brunius C, Shi L and Landberg R. Within and between batch correction of LC-MS metabolomics data. Submitted manuscript.*

## Version history
version | date  | comment
:-----: | :---: | :------
0.1.7 | 16-04-12 | Added `data` statement to `workflow.R` and added population normalization to `refCorr`.
0.1.6 | 16-02-03 | Added `batchComb` function to extract features present in multiple batches. Updated `workflow.R` <- Removed old data-raw (included in `workflow.R`). Moved from GitHub to GitLab
0.1.5 | 15-12-17 | Updated `workflow.R`. Added `grabWrapBA` for grabbing using batch-aligned peaktable. Added batchdata. Updated peaktable data.
0.1.4 | 15-12-03 | Added `workflow.R` under `inst/`, roxygenised data and updated README.
0.1.3 | 15-12-03 | Added data and data-raw
0.1.2 | | `batchFlag()`: externalised peakInfo (updated `peakInfo()` to include starting character) 
0.1.1 | | Roxygenised. Successful build.
0.1.0 | | Functions in place. Not fully roxygenised.

