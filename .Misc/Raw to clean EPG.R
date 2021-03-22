# Import raw treatment data -----------------------------------------------
library(data.table)
EPG.raw<- data.table(read.csv(".Misc/Parasite.csv"))

# rename consistently season names
EPG.raw[season=="2.win",season:="2.winter"]

# clean the season factor's levels
# droplevels(EPG.raw$season) # not needed since data.table update
EPG.raw$date<- as.Date(EPG.raw$date,origin="1970-01-01")

# consistently format chronological season's name
EPG.raw[,seasonchron:=paste0(year(date),"-",season)]

# subset raw columns into those useful for FECR script
ordered.col.to.keep<- c("date","season","seasonchron","sample","id","tvnt","treatment","species","epg")
EPG.clean<- EPG.raw[,ordered.col.to.keep,with=FALSE][order(seasonchron)]

# make sample date will be handled as a date afterward, and year as explicit data
EPG.clean$date<- as.Date(EPG.clean$date,origin="1970-01-01")
EPG.clean$year<- year(EPG.clean$date)

# format it into regular data.frame
EPG<- data.frame(EPG.clean)

# export into ./data/
usethis::use_data(EPG,overwrite = TRUE)
