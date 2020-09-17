# Swing Ratio Analysis

Contains data and analysis code for a study on swing ratios in jazz solo improvisation.

## Usage
The main analysis file is `swing_ratio_anaylsis.R`. Source file (probably, you need to install some missing R packages first), and run `setup_workspace(recalc_GMM = T)`. The main data object is `srw_f`, Gaussian mixture models for a variety of subgroups can be found in `solo_gmm`, `performer_gmm`, `swing_gmm`, `tempo_gmm`, `decade_gmm`.  

Raw data is in `data/swing_ratios_raw.csv`. Cleaned, filtered, and annotated data is in `data/swing_ratios_cleaned_with_loudness.csv`. This files contains only events from swing triples, and has also loudness and pitch as well as some useful metadata annotations.

## Citation


## Acknowledgments

 
## Implementation notes
This repo is still under development. Please check again for updates.
