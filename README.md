# Genetic overlap between RA and MI

This repository hosts scripts and documents, all related to the study described in [Sysojev et al. (2024)](https://acrjournals.onlinelibrary.wiley.com/doi/10.1002/art.42918). Worth noting, is that it does not contain any type of raw data used within the project, whatsoever.

# 1. SCRIPTS

The project is divided into four primary scripts, with supporting `misc` scripts underlying them. Script `1*.sh` performs the pre-processing of the fully imputed genotype data by performing quality control in each of the underlying cohorts and merging of the post-QC data into a singular set of genotype data. This data is then taken forward into script `2*.sh`, wherein the genome-wide association analysis is conducted, and GWAS summary statistic data is obtained. In script `3*.sh`, the GWAS data is combined with publicly available GWAS data on phenotypes of intrest, in which genome-wide genetic correlation is estimated through LDSC.

Script `4*.R` contains code for performing analyses with LAVA, for the estimation of local genetic correlations. Here, there may be inconsistencies in terms of filepaths as it was built on my Windows system, whereas the previous scripts were built on the Linux server. This was done to exploit the greater set of available cores, to improve speed in parallelization. While this script *works*, I do recommend exerting care if aiming to look to it for inspiration for similar analyses, as it is written in a somewhat ad-hoc way.

A few of the `misc` scripts have no super script which controls them: these are primarily scripts that extract data for tables employed within the paper. As such, they can be considered part of a non-existent `5*.R` script that pulls data from the registers and processes it for descriptive statistics.

# 2. DOCS

Various types of documents and documentation related to the project, most likely being self-explanatory. The `misc` folder contains a few minor files text files with notes relating to issues and challenges within the project. While not key, I've kept them as it is possible that others may have use of them. Other folders have been cleaned and made to only include the latest version of files, to avoid clutter and unnecessary repetition - earlier versions can be made available upon request, of course.