# Import datasets (epg and treatment data) --------------------------------
library(data.table)
EPG<- data.table(read.csv("data/Parasite.csv"))
Treatment<- data.table(read.csv("data/Treatment.csv"))

EPG$date<- as.Date(EPG$date,origin="1970-01-01")
EPG$year<- year(EPG$date)
Treatment$Trtdate<- as.Date(Treatment$Trtdate,origin="1970-01-01")

closest.date<- function(D,D.list){
  D.list[which.min(D,D.list)]
}


there.is.date<- sapply(seq_len(nrow(EPG)),
                       function(i) {
                         year.epg<- EPG$year[i];season.epg<- EPG$season[i];
                         id.epg<- EPG$id[i];treatment<- EPG$treatment[i];
                         date<- EPG$date[i];
                         
                         date.list<- Treatment[id==id.epg & year==year.epg & season==season.epg]$Trtdate
                         cat(paste0("\rcalculating date ",i,"/",nrow(EPG)))
                         closest<- closest.date(date = date,date.list = date.list)
                         length(closest)!=0
                       }
)

EPG<- EPG[there.is.date]

EPG$Trtdate<- as.Date(
  sapply(seq_len(nrow(EPG)),
         function(i) {
           year.epg<- EPG$year[i];season.epg<- EPG$season[i];
           id.epg<- EPG$id[i];treatment<- EPG$treatment[i];
           D<- EPG$date[i];
           
           D.list<- max(as.Date(Treatment[id==id.epg & year==year.epg & season==season.epg]$Trtdate))
           cat(paste0("\rcalculating date ",i,"/",nrow(EPG)))
           closest.date(D,D.list)
         }
  ),
  origin="1970-01-01"
)

EPG<- EPG[EPG$date>=(EPG$Trtdate-3*7) & EPG$date<=(EPG$Trtdate+3*7)]   # confirm timing for inclusion

EPG$period<- ifelse(EPG$date<EPG$Trtdate,"before","after")   # what about on the day => "before" to confirm.



