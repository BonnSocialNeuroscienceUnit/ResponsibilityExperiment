# README file for Responsibility study

## Contents
This repository contains Matlab code (.m) and data files (.mat, .csv, .nii) required to reproduce the results of the study entitled "A neural mechanism of social responsibility", by Maria Gaedeke, Tom Willems, Omar Salah Ahmed, Bernd Weber, Rene Hurlemann, and Johannes Schultz. A preprint of this study can be found on bioRxiv (Link: __https://doi.org/10.1101/2020.05.25.107300__). The fMRI data of the study are available on OpenNeuro (Link: __https://openneuro.org/datasets/ds005588__).

### Structure of this repository
__BehaviouralData__ contains the raw behavioural data (.mat format) from the two studies reported in the manuscript.
__Code__ contains Matlab code to generate the figures and reproduce the analyses (see below for details). Within Code, the folders Figures and csv contain outputs from the analyses, fMRI contains fMRI-specific Matlab scripts, and bin contains helper functions.
__fMRIresults__ contains Nifti files resulting from the analysis of the fMRI data and used in the Figures.

## How to replicate the behavioural analyses and figures
To replicate the results Figures using Matlab:
First, download the repository and add all folders to the path.
Then run:

- master\_behavAnalysis.m to generate Figures 2 and 3, and Supplementary Figure 1
- master\_fMRI\_showResults to generate Figures 4 and 5, and Supplementary Figure 2

## How to replicate the fMRI analysis
To replicate the fMRI results using SPM12 in Matlab based on the raw data, download the raw data from OpenNeuro (Link: __https://openneuro.org/datasets/ds005588__). Then, use the following scripts located in Code/fMRI:

- a\_preprocessing.m for realignment / motion correct, normalization/unwarping, and smoothing
- b\_FFX.m to fit individual general linear models to the BOLD data
- c\_gPPI.m to run functional connectivity analyses using the gPPI toolbox (download on https://www.nitrc.org/projects/gppi)
- d\_RFX\_fullFact.m to run most group-level "random effects" models
- e\_RFX\_pairedTtest\_guiltEffect.m for a specific group-level test
- f\_apply\_YuKoban\_guilt\_signature\_map.m to compare our findings to previously published results (uses code from https://github.com/canlab/Neuroimaging\_Pattern\_Masks)

## Contact
johannes.schultz@ukbonn.de
