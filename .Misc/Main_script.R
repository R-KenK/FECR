# Import datasets (epg and treatment data) --------------------------------
source("R/Misc_tools.R")
source("R/FECR.R")

load(file = "data/EPG.rda")
load(file = "data/Treatment.rda")

EPG$date<- as.Date(EPG$date,origin="1970-01-01")
EPG$year<- as.numeric(substr(EPG$date,1,4))
Treatment$Trtdate<- as.Date(Treatment$Trtdate,origin="1970-01-01")

there.is.date<- sapply(seq_len(nrow(EPG)),
                       function(i) {
                         year.epg<- EPG$year[i];season.epg<- as.character(EPG$season[i]);
                         id.epg<- as.character(EPG$id[i]);treatment<- EPG$treatment[i];
                         D<- EPG$date[i];

                         D.list<- Treatment[Treatment$id==id.epg & Treatment$year==year.epg & Treatment$season==season.epg,]$Trtdate
                         cat(paste0("\rcalculating date ",i,"/",nrow(EPG)))
                         length(closest.date(D,D.list))!=0
                       }
)

EPG<- EPG[there.is.date,]

EPG$Trtdate<- as.Date(
  sapply(seq_len(nrow(EPG)),
         function(i) {
           year.epg<- EPG$year[i];season.epg<- as.character(EPG$season[i]);
           id.epg<- as.character(EPG$id[i]);treatment<- EPG$treatment[i];
           D<- EPG$date[i];

           D.list<- max(as.Date(Treatment[Treatment$id==id.epg & Treatment$year==year.epg & Treatment$season==season.epg,]$Trtdate))
           cat(paste0("\rcalculating date ",i,"/",nrow(EPG)))
           closest.date(D,D.list)
         }
  ),
  origin="1970-01-01"
)

EPG<- EPG[EPG$date>=(EPG$Trtdate-3*7) & EPG$date<=(EPG$Trtdate+3*7),]   # confirm timing for inclusion

EPG$period<- ifelse(EPG$date<EPG$Trtdate,"before","after")   # what about on the day => "before" to confirm.



