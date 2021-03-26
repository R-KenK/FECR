# Import required functions and files --------------------------------
source("R/seasonchron_tools.R")

# Import datasets (epg and treatment datasets from .csv, exported into the FECR package) --------------------------------
load(file = "data/EPG.rda")
setDT(EPG)

# match samples to the season(chron) defined in seasonchron.lookup
seasonchron.lookup <- create_seaonchron.lookup()
EPG <- match_epg_seasonchron(epg.df = EPG,
                         seasonchron.lookup = change_seasonchron_start.end(seasonchron.lookup,"2014-2.spring","2014-01-15")
                         # seasonchron.lookup = seasonchron.lookup
)
