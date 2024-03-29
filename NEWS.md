# FECR 5.0.0
* Added a new FECR calculation method, "MacIntosh2", compatible with multiple epg value per
individual.
* Method "MacIntosh" has therefore been renamed "MacIntosh1"
* Non-parametric, the method repeats many times the followings:    
    * randomly sample one epg value per individual
    * on the randomly sampled epg values, previously introduced method "MacIntosh1" is applied.
    * the mean FECR value is reported, as well as the confidence interval of choice.

# FECR 4.0.1
* Added `"pb"` argument to `FECR()` to not display a progress bar

# FECR 4.0.0
* cleaned package's structure, preparing for clean release
* moved check and averaging functions to separate files, moves conditional averaging to a wrapper function

# FECR 3.0.1
* updated NEWS.md

# FECR 3.0.0
* added confidence interval capabilities to FECR4-6.
* implemented CIs in FECR calculation and FECR barplot.

# FECR 2.1.1
* removed left overs referring to now unused variable (fixed Andrew's error).

# FECR 2.1.0
* added few epg modelling stats and few exploratory plots.

# FECR 2.0.0
* revamped data handling with explicit raw to clean dataset scripts.
* this allowed the addition of a new season.
* script homogenization and revamping.

# FECR 1.0.0
* Added a `NEWS.md` file to track changes to the package.
* First attempt at data handling and script backbone to compute FECRs.
* Include base FECR() and functions in Misc_tools.
