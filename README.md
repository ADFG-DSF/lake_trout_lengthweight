# Length-Weight Relationship and Asymptotic Length for Lake Trout

The purpose of this modeling exercise was to provide Alaska-specific estimates of the length-weight relationship in lake trout, as well as the possible relationship between lake area and asymptotic weight, as a means to further refine the modeling framework described in Lester et al. (yyyy) for estimation of Maximum Sustained Yield of lake trout in Alaska.

Data were provided from [source??] and compiled by Region III biologists as relevant.  Data represent a total of:
* 82 lakes
* 137 sampling events from 1960-2022
* 34,362 fish sampled
* 34,107 length measurements
* 3,894 weight measurements
* 3,887 paired length x weight measurements from
  - 24 lakes
  - 57 sampling events from 1960-2022

## Repository contents

### flat_data/

* `lake_morphometry.csv` gives lake-level data (area, latitude, etc)

* `length_weight.csv` gives fish-level data (length, weight, etc)

### R/

* `1_laketrout_lwdata.R` reads data and does some basic cleaning for the purpose of length-weight relationship

* `2_laketrout_lwmodels.R` fits an exhaustive set of length-weight models (log-log regression), with all possible combinations of slope and intercept parameters:
  - common for all lakes
  - separate for all lakes but hierarchically distributed
  - trending with log-area
  - trending with latitude

## Reference

Lester, et al. citation
