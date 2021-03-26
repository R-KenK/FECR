# Import required functions and files --------------------------------
source("R/seasonchron_tools.R")

# Import datasets (epg and treatment datasets from .csv, exported into the FECR package) --------------------------------
load(file = "data/Treatment.rda")
Treatment
setDT(Treatment)

# match samples to the season(chron) defined in seasonchron.lookup
seasonchron.lookup <- create_seaonchron.lookup()
Treatment <- match_treatment_seasonchron(Treatment.df = Treatment,
                         seasonchron.lookup = change_seasonchron_start.end(seasonchron.lookup,"2014-2.spring","2014-01-15")
                         # seasonchron.lookup = seasonchron.lookup
)
