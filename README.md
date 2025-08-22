# Estimating Asymptotic Weight for Lake Trout

Application of the Lester Model to lake trout in Alaska lakes may be substantially improved by direct estimation of asymptotic weight W_∞ when possible. In the computational framework presented by Lester, et al., W_∞ is estimated by means of a sequence of relationships: first, the relationship observed between lake area and asymptotic length observed in Canada; then, the aggregate relationship between length and weight observed in Canada.

We have data at a finer resolution: for 41 lakes, we have sampling data on 30,000+ individual lake trout (4,772 weights; 32,574 lengths; 1,517 ages), which can be used to gain inference on pairwise relationships between variables (4,764 paired length x weight observations from 26 lakes; 1,517 paired length x age observations from 19 lakes). Not only can the relationships between length and weight and between lake area and asymptotic length be estimated for Alaska, they can be estimated specifically for each of the 41 lakes and applied to additional lakes as appropriate, and parameters reported by Lester, et. al. may be incorporated as priors.

Implementation of an integrated model as outlined below allows simultaneous use of all available information, as well as appropriate propagation of all uncertainty in estimation. This should represent an improvement over a previous modeling exercise that estimated each relationship in sequence and then employed heuristic rules to select among available methods for estimating asymptotic size: rather, all available information is used simultaneously.


## Repository contents

### flat_data/

* `lake_morphometry3.csv` gives lake-level data (area, latitude, etc)

* `length_weight4.csv` gives fish-level data (length, weight, etc)

* all other .csv files are no longer current and will be cleaned out soon.

### R/

* `5_integrated_Winf.R` is a standalone script that fits the integrated Bayesian model

* all other .R scripts in this folder are analyses that are no longer current!

### integrated model equations.Rmd 

This Markdown document outlines and describes the current version of the integrated Bayesian model

### Winf_estimates.csv

This spreadsheet summarizes model inferences and data inputs for all lakes

## Reference

Lester, et al. citation
