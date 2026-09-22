# Using Gaussian processes to quantify fine-scale heterogeneity in wildland fuels

## Abstract

A suite of hierarchical spatial statistical models using Gaussian process spatial 
random effects is proposed to quantify fine-scale spatial patterns in the fuel 
load of several wildland fuel components (duff, litter, and fine woody debris). 
A sampling protocol and associated statistical models are described and applied 
in a Sierra Nevada mixed conifer forest affected by extensive drought caused 
tree mortality. Model validation reveals varying performance in three tasks: 
1) Making pointwise predictions, 2) reproducing the distribution of fine-scale 
load observations, and 3) reproducing the distribution of coarse-scale mean 
fuel loads. Models for the depth of duff, depth of litter, and count of fine 
woody debris particles perform well in all three tasks and parameter estimates 
are well informed by the data. The sampling protocol and statistical analysis 
described in this study enable quantitative description and simulation of the 
fine-scale spatial patterns of fuel loads, an important step towards the 
creation of fuel model inputs for next-generation fire spread modeling with 
fine-scale spatial processes. 


## This repo

- **00-R/** contains generic functions for using Brown's transects data to estimate 
fuel loads
- **01-context/** contains some background and planning documents
- **02-data/** contains the raw data and intermediate data products fed to the 
stan models. 
  - The actual model results are deposited in 02-data/06-results, but 
are too large to commit to git. See notes below to reproduce them. 
  - There are also a few large geospatial datasets excluded from the commit of 
  `02-data/00-source/`. These are described in the manuscript and available from 
  their public sources.
- **03-analysis/** contains the core analysis scripts for ingesting the data, 
defining the models, estimating the parameters, and analyzing the results.
  - R scripts are numbered according to their order in the data pipeline
  - xx- and zz- files are work branches that did not wind up in the final paper
  - Stan models (both stan code and compiled executables) are labeled accordingly. 
  See `03-analysis/05-estimate_parameters.R` for code pointing to the specific 
  stan code used for the paper.
- **04-communication/** contains figures, tables, and the manuscript.

## Reproducing fitted model results

The stan model objects and draws .csvs are produced by `03-analysis/08-estimate_parameters.R`, 
which loads the data objects and runs the stan models to estimate parameters. 
These stan model objects and the draws .csvs are too large to commit to github.
Scripts downstream of `03-analysis/08-estimate_parameters.R` require the stan 
model objects. To reproduce the figures and tables included in the manuscript, 
you will need to re-run `03-analysis/08-estimate_parameters.R` to produce those 
files locally. 
