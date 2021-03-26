# Import raw treatment data -----------------------------------------------
library(data.table)
Trt.raw<- data.table(read.csv(".Misc/Raw Treatment.csv"))

# rename consistently season names
Trt.raw[season=="fall",season:="4.fall"]
Trt.raw[season=="winter",season:="1.winter"]
Trt.raw[season=="spring",season:="2.spring"]
Trt.raw[season=="summer",season:="3.summer"]

# clean the season factor's levels
# droplevels(Trt.raw$season) # not needed since data.table update
Trt.raw$date<- round(as.Date(Trt.raw$date,origin="1970-01-01"))

# consistently format chronological season's name
Trt.raw[,seasonchron:=paste0(year,"-",season)]

# subset raw columns into those useful for FECR script
ordered.col.to.keep<- c("year","season","seasonchron","id","drug","dose_no","date")
Trt.clean<- Trt.raw[,ordered.col.to.keep,with=FALSE][order(seasonchron)]

# determine controls' average treatment date (of the relevant seasonchron) for each dose number
Trt.control<- Trt.raw[,.(date=round(mean(date)),id="control",drug="control"),by=.(year,season,seasonchron,dose_no)][,ordered.col.to.keep,with=FALSE]

# add controls' average date to the treatment data of the treated group
Trt.clean<- rbind(Trt.clean,Trt.control)[order(seasonchron)]

# rename the date column (easier to know which date it is about in the script), and make sure it'll be handled as a date afterward
colnames(Trt.clean)[7]<- "Trtdate"
Treatment$Trtdate<- as.Date(Treatment$Trtdate,origin="1970-01-01")

# format it into regular data.frame
Treatment<- data.frame(Trt.clean)

# export into ./data/
usethis::use_data(Treatment,overwrite = TRUE)
