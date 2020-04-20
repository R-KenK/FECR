# Import required functions --------------------------------
source("R/Misc_tools.R")
source("R/FECR.R")

# Import datasets (epg and treatment datasets from .csv, exported into the FECR package) --------------------------------
load(file = "data/EPG.rda")
load(file = "data/Treatment.rda")

# Keep only oes and trich
EPG<- EPG[EPG$species %in% c("oes","trich"),]
EPG$species<- droplevels(EPG$species)

# Identify and subset epg data with treatment for the specific season
id.or.control<- ifelse(EPG$tvnt=="treated",as.character(EPG$id),"control")
there.is.date<- sapply(seq_len(nrow(EPG)),
                       function(i) {
                         seasonchron.epg<- as.character(EPG$seasonchron[i]);
                         id.epg<- as.character(id.or.control[i]);
                         treatment<- EPG$treatment[i];D<- EPG$date[i];

                         D.list<- Treatment[Treatment$id==id.epg & Treatment$seasonchron==seasonchron.epg,]$Trtdate
                         cat(paste0("\rcalculating date ",i,"/",nrow(EPG)))
                         length(D.list)!=0
                       }
)
EPG<- EPG[there.is.date,]

# Match epg data with the relevant individual or average treatment date
EPG$Trtdate<- as.Date(
  sapply(seq_len(nrow(EPG)),
         function(i) {
           seasonchron.epg<- as.character(EPG$seasonchron[i]);
           id.epg<- as.character(EPG$id[i]);tvnt.epg<- as.character(EPG$tvnt[i]);
           treatment<- EPG$treatment[i];D<- EPG$date[i];

           D.list<- max(as.Date(Treatment[Treatment$id==ifelse(tvnt.epg=="treated",id.epg,"control") & Treatment$seasonchron==seasonchron.epg,]$Trtdate))
           cat(paste0("\rcalculating date ",i,"/",nrow(EPG)))
           closest.date(D,D.list)
         }
  ),
  origin="1970-01-01"
)

# Keep only epg data when the collection date is in the relevant time window with respect to treatment date
EPG.all<- EPG
EPG.3w<- EPG[EPG$date>=(EPG$Trtdate-3*7) & EPG$date<=(EPG$Trtdate+3*7),]   # confirm timing for inclusion

# Determine before and after periods based on matched treatment date (revamping the previous sort-of-manual period attribution)
EPG.3w$period<- ifelse(EPG.3w$date<EPG.3w$Trtdate,"before","after")   # what about on the day => "before" to confirm.

# Subsetting columns of interest
FECR.df<- EPG.3w[,c("date","season","seasonchron","sample","id","tvnt","Trtdate","period","species","epg")]

# Subset epg data with only samples from individual that were sampled before and after treatment
FECR.df<-  keep_with_before.after(FECR.df) # NOT SURE ABOUT INCLUDING THIS LINE OR NOT: should we remove individuals without data in both before and after, when considering mean values?

# Calculate mean (and sd) epg values of the different C1,C2,T1,T2 of each season
FECR.group<- as.data.frame(data.table::data.table(FECR.df)[,.(date=mean(date),epg=mean(epg),sd.epg=sd(epg)),by=.(seasonchron,tvnt,season,period,species)][order(seasonchron)])

# Same, but at the individual level
FECR.ind<- as.data.frame(data.table::data.table(FECR.df)[,.(date=mean(date),epg=mean(epg),sd.epg=sd(epg)),by=.(seasonchron,tvnt,id,season,period,species)][!(tvnt=="treated"&period=="before"&epg==0)][!(tvnt=="control"&period=="after"&epg==0)])
FECR.ind<-  keep_with_before.after(FECR.ind)

FECR.comp<- rbind_lapply(unique(FECR.df$species),
                         function(sp){
                           FECR.group<- data.table::data.table(FECR.group);FECR.ind<- data.table::data.table(FECR.ind);
                           rbind_lapply(unique(FECR.df$seasonchron),
                                        function(s){
                                          C1.mean<- FECR.group[seasonchron==s & tvnt=="control" & period=="before" & species==sp]$epg
                                          C2.mean<- FECR.group[seasonchron==s & tvnt=="control" & period=="after" & species==sp]$epg
                                          T1.mean<- FECR.group[seasonchron==s & tvnt=="treated" & period=="before" & species==sp]$epg
                                          T2.mean<- FECR.group[seasonchron==s & tvnt=="treated" & period=="after" & species==sp]$epg
                                          
                                          if(any(length(C1.mean)==0,length(C2.mean)==0,length(T1.mean)==0,length(T2.mean)==0)) {return(NULL)}
                                          
                                          C1.ind.list<- C2.ind.list<- rand.ind.list(FECR.ind,s = s,group = "control",per = "before",sp = sp)
                                          T1.ind.list<- T2.ind.list<- rand.ind.list(FECR.ind,s = s,group = "treated",per = "before",sp = sp)
                                          
                                          C1.sub<- FECR.ind[seasonchron==s&period == "before"&species==sp][id %in% C1.ind.list]$epg
                                          C2.sub<- FECR.ind[seasonchron==s&period == "after"&species==sp][id %in% C1.ind.list]$epg
                                          T1.sub<- FECR.ind[seasonchron==s&period == "before"&species==sp][id %in% T1.ind.list]$epg
                                          T2.sub<- FECR.ind[seasonchron==s&period == "after"&species==sp][id %in% T1.ind.list]$epg
                                          
                                          C1<- FECR.ind[seasonchron==s&tvnt=="control"&period == "before"&species==sp]$epg
                                          C2<- FECR.ind[seasonchron==s&tvnt=="control"&period == "after"&species==sp]$epg
                                          T1<- FECR.ind[seasonchron==s&tvnt=="treated"&period == "before"&species==sp]$epg
                                          T2<- FECR.ind[seasonchron==s&tvnt=="treated"&period == "after"&species==sp]$epg
                                          
                                          data.frame(species=sp,seasonchron=s,
                                                     FECR1=FECR(T1 = T1.mean,T2 = T2.mean,method = "Kochapakdee"),
                                                     FECR2=FECR(T1 = T1.mean,T2 = T2.mean,C1 = C1.mean,C2 = C2.mean,method = "Dash"),
                                                     FECR3=FECR(T2 = T2.mean,C2 = C2.mean,method = "Coles"),
                                                     FECR4=FECR(T1 = T1.sub,T2 = T2.sub,method = "Cabaret1"),
                                                     FECR5=FECR(T1 = T1.sub,T2 = T2.sub,C1 = C1.sub,C2 = C2.sub,method = "Cabaret2"),
                                                     FECR6=FECR(T1 = T1,T2 = T2,C1 = C1,C2 = C2,method = "MacIntosh")
                                          )
                                        }
                           )
                         }
)
# Display the table of FECR values
FECR.comp


library(ggplot2)
# raw plotting of the data
EPG.all$treatment<- factor(EPG.all$treatment,levels = c("before","after"),ordered = TRUE)
ggplot(EPG.all,aes(treatment,epg,colour=tvnt))+
  facet_wrap(.~seasonchron+species+tvnt,ncol=4,scales = "free_y")+
  geom_jitter(alpha=0.3)+geom_boxplot(alpha=0.75)+
  theme_bw()

# Effect of including only +/- 3 weeks or raw data
EPG.3w$treatment<- factor(EPG.3w$treatment,levels = c("before","after"),ordered = TRUE)
ggplot(EPG.3w,aes(treatment,epg,colour=tvnt))+
  facet_wrap(.~seasonchron+species+tvnt,ncol=4,scales = "free_y")+
  geom_jitter(alpha=0.3)+geom_boxplot(alpha=0.75)+
  theme_bw()

# Group epgs
FECR.group$period<- factor(FECR.group$period,levels = c("before","after"),ordered = TRUE)
ggplot(FECR.group[FECR.group$seasonchron!="2014-2.winter",],aes(period,epg,fill=period))+
  facet_wrap(seasonchron~species+tvnt,ncol=4,scales = "free_y")+
  geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin=epg,ymax=epg+sd.epg),width=0.25)+
  geom_bar(stat = "identity")+
  theme_bw()

# individual epgs (note that no notion of intra-individual variation appears within C1,C2,T1,T2)
FECR.ind$period<- factor(FECR.ind$period,levels = c("before","after"),ordered = TRUE)
ggplot(FECR.ind[FECR.ind$seasonchron!="2014-2.winter",],aes(period,epg,colour=period,group=id))+
  facet_wrap(seasonchron~species+tvnt,ncol=4,scales = "free_y")+
  geom_hline(yintercept = 0)+
  geom_line(alpha=0.5)+
  geom_point()+
  theme_bw()

# FECRs plots
FECR.long<- reshape2::melt(FECR.comp,id.vars = c("species","seasonchron"),variable.name = "method",value.name = "FECR")
ggplot(FECR.long,aes(method,FECR,fill=method))+
  facet_wrap(seasonchron~species,ncol=2,scales = "free_y")+
  geom_hline(yintercept = 0)+geom_hline(yintercept = 100,lty="dashed",colour="grey50")+
  geom_bar(stat = "identity")+
  theme_bw()

# quick attempt at having well-fitted models, with lowest AIC by tinkering roughly...
library(glmmTMB)
FECR.oes.glmm<- glmmTMB(round(epg)~period+period:tvnt+seasonchron:period+(1|id),
                        ziformula = ~period+period:tvnt+seasonchron,
                        dispformula = ~1,
                        data = FECR.df[FECR.df$species=="oes"&FECR.df$seasonchron!="2014-2.winter",],family = "nbinom2")
# model fit looks ok, and this combination of fixed/random/zi/disp effects yielded the lowest AICs
plot(DHARMa::simulateResiduals(FECR.oes.glmm))
summary(FECR.oes.glmm)
FECR.oes.null<- glmmTMB(round(epg)~1+(1|id),
                        data = FECR.df[FECR.df$species=="oes"&FECR.df$seasonchron!="2014-2.winter",],family = "nbinom2")

# this would have been interesting...
summary(FECR.oes.null)
MuMIn::r.squaredGLMM(FECR.oes.glmm,FECR.oes.null) # ... but this is unfortunately not available :(

FECR.trich.glmm<- glmmTMB(round(epg)~period+period:tvnt+seasonchron:period+(1|id),
                          ziformula = ~period+period:tvnt+seasonchron,
                          dispformula = ~1,
                          FECR.df[FECR.df$species=="trich"&FECR.df$seasonchron!="2014-2.winter",],family = "nbinom2")
# model fit looks ok, and this combination of fixed/random/zi/disp effects yielded the lowest AICs
plot(DHARMa::simulateResiduals(FECR.trich.glmm))
summary(FECR.trich.glmm)
FECR.trich.null<- glmmTMB(round(epg)~1+(1|id),
                          FECR.df[FECR.df$species=="trich"&FECR.df$seasonchron!="2014-2.winter",],family = "nbinom2")
# this would have been interesting...
summary(FECR.trich.null)
MuMIn::r.squaredGLMM(FECR.trich.glmm,FECR.trich.null) # ... but this is unfortunately not available :(

# but for which the across season variability impact results (and model predictions) heavily (careful, wrongly ordered before and after...)
plot(ggeffects::ggpredict(FECR.oes.glmm,c("tvnt","period","seasonchron")))
plot(ggeffects::ggpredict(FECR.trich.glmm,c("tvnt","period","seasonchron")))
