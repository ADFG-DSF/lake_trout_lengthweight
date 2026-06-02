# Estimating Asymptotic Weight for Lake Trout

Application of the Lester Model to lake trout in Alaska lakes may be 
substantially improved by direct estimation of asymptotic weight 
$W_\infty$ when possible. In the computational framework presented by 
Lester, et al.,$W_\infty$ is estimated by means of a sequence of 
relationships: first, the relationship observed between lake area 
and asymptotic length observed in Canada; then, the aggregate 
relationship between length and weight observed in Canada.

We have data at a finer resolution: for 44 lakes, we have sampling 
data on 40,000+ individual lake trout (5,695 weights; 36,041 lengths; 
1,494 ages), which can be used to gain inference on pairwise 
relationships between variables (5,695 paired length x weight 
observations from 29 lakes; 1,494 paired length x age observations 
from 19 lakes). Not only can the relationships between length and 
weight and between lake area and asymptotic length be estimated for 
Alaska, they can be estimated specifically for each of the 44 lakes 
and applied to additional lakes as appropriate, and parameters 
reported by Lester, et. al. may be incorporated as priors.

Implementation of an integrated model as outlined below allows 
simultaneous use of all available information, as well as 
appropriate propagation of all uncertainty in estimation. 


## Repository contents

### /flat_data/

* `lake_morphometry_25_12_23.csv` gives lake-level data (area, latitude, etc)

* `length_weight_25_12_23.csv` gives fish-level data (length, weight, etc)

* all other .csv files are no longer current but are retained for checks.

### /R/

The R scripts are intended to be run more or less in sequence as appropriate.

* `0_data_update_checks.R` runs a sequence of checks that can be run whenever
the lake-level or fish-level datasets are updated, hopefully alerting the 
user to any weirdness.

* `1_laketrout_data_import.R` performs all data import AND the current 
data filters to create the dataset that is assumed to be clean and will
be used going forward.  Outputs of this script will be loaded by latter
scripts in the sequence.

* `2_integrated_model_run.R` is specific to running the integrated Bayesian
model.  It bundles the input data, defines the model itself, runs the model,
saves the output, and that's it.

* `2a_LWmodel_exploration.R` performs a sequence of checks to explore several
candidate versions of the weight~length component of the model.  It got
a little out of control.

* `3_plotting_model_run.R` loads the output created by script 2 and plots
the output every way I could think of.  It got even more out of control.

* **/R/further_model_exploration/** is some stuff that further supplements
script 2a, and will hopefully be cleaned out a bit.

### /Rdata/

This is where .Rdata files to be loaded later will live

### Report_draft.Rmd

This is the current working draft of the overall report!  The output Word
file will need lots of formatting and the contents will most likely be
pasted into a ADFG report template.

### flowchart1.pptx and flowchart1.jpg

The model in conceptual flowchart form.  This should probably not live in the root directory

### Appendix_LakePlots.Rmd

This creates many summary plots for all lakes of interest, and is 
kept contained to its own file, perhaps as a warning to others.

### Winf_estimates.csv

This spreadsheet summarizes model inferences and data inputs for all lakes.

## Reference

Lester, et al. citation
