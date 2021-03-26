
# data exploratory & visualization ----------------------------------------

FECR.df
EPG.all

library(ggplot2)

EPG.all$treatment <- factor(EPG.all$treatment,levels = c("before","after"))
ggplot(EPG.all,aes(interaction(treatment,tvnt),epg,fill = treatment,colour = treatment))+
  facet_wrap(species~seasonchron,scales = "free",nrow = 2)+
  geom_line(aes(group = id),colour = 'grey80')+
  geom_point(alpha=.3)+
  geom_boxplot(colour = "black",alpha=.3)+
  theme_bw()

EPG.all$period<- ifelse(EPG.all$date < EPG.all$Trtdate,"before","after")

EPG.all$treatment == EPG.all$period


min.max.raw.trt <- Trt.raw[,.(min = min(date),max = max(date)),by = .(year,season)]
min.max.clean.trt <- Treatment[,.(min = min(Trtdate),max = max(Trtdate)),by = .(year,season)]

min.max.raw.epg <- EPG.raw[,.(min = min(date),max = max(date)),by = .(year(date),season)]
min.max.clean.epg <- EPG[,.(min = min(date),max = max(date)),by = .(year,season)]


date.check <- rbind(
  # data.table(type = "trt",processed = "raw",min.max.raw.trt),
  data.table(type = "trt",processed = "clean",min.max.clean.trt),
  # data.table(type = "epg",processed = "raw",min.max.raw.epg),
  data.table(type = "epg",processed = "clean",min.max.clean.epg)
)

date.check.long <- melt.data.table(date.check,measure.vars = c("min","max"),value.name = "date")

ggplot(date.check.long,aes(date,type,colour= variable))+
  facet_grid(year+season~processed)+
  geom_jitter(data = EPG,aes(date,y = "epg"),colour = "grey50",alpha = .05,height = .2)+
  geom_vline(aes(xintercept = date,colour= variable),lty = "dashed" )+
  geom_vline(data = date.check.long[type == "trt",.(date = mean(date)),by = .(year,season)],aes(xintercept = date),colour= "red")+
  scale_x_date(date_breaks = "1 month")+
  geom_point()+
  theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.minor = element_blank())


  # main script redux draft -------------------------------------------------

# Identify and subset epg data with treatment for the specific season
id.or.control<- ifelse(EPG$tvnt=="treated",EPG$id,"control")
there.is.date<- sapply(seq_len(nrow(EPG)),
                       function(i) {
                         seasonchron.epg<- EPG$seasonchron[i];
                         id.epg<- id.or.control[i];
                         treatment<- EPG$treatment[i];D<- EPG$date[i];
                         
                         D.list<- Treatment[Treatment$id==id.epg & Treatment$seasonchron==seasonchron.epg,]$Trtdate
                         cat(paste0("\rcalculating date ",i,"/",nrow(EPG)))
                         length(D.list)!=0
                       }
)
EPG<- EPG[there.is.date,]



# functions ---------------------------------------------------------------


# Import required functions -----------------------------------------------

library(ggplot2)
library(data.table)
source("R/Misc_tools.R")
source("R/seasonchron_tools.R")
source("R/subsetting_tools.R")
source("R/FECR.R")
source("R/FECR_tools.R")

# Import datasets (epg and treatment datasets from .csv, exported into the FECR package) --------------------------------

load(file = "data/EPG.rda")
load(file = "data/Treatment.rda")

# visualization of the timings and data inclusion --------------------

## raw data
ggplot(
  # melt(change_seasonchron_start.end(seasonchron.to.change = "2014-2.spring",new.start = "2014-01-15"),
  melt(create_seaonchron.lookup(),
       measure.vars = c("start","end"),value.name = "date"),
  aes(date,seasonchron,colour = variable))+
  geom_point(shape = 22,aes(fill = variable),colour = "grey20", size = 2)+
  geom_jitter(
    data = EPG[EPG$species == "oes",],
    height = .25,aes(date,seasonchron,colour = treatment),alpha = .3)+
  geom_jitter(
    data = Treatment,
    aes(Trtdate,seasonchron),colour = "grey30",fill = "tomato",shape = 21,size = 2,alpha = .2,height = .1)+
  scale_x_date(date_breaks = "1 month")+
  theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.minor = element_blank())


# matching epg and treatment data -----------------------------------------
merged.data <- merge_epg_and_treatment_data(EPG,Treatment,
                             # seasonchron.to.change = "2014-2.spring",new.start = "2014-01-15"
                             seasonchron.to.change = NULL
)
EPG <- merged.data$epg.df
Treatment <- merged.data$trt.df

## no change in seasonchron data
ggplot(
  # melt(change_seasonchron_start.end(seasonchron.to.change = "2014-2.spring",new.start = "2014-01-15"),
  melt(create_seaonchron.lookup(),
       measure.vars = c("start","end"),value.name = "date"),
  aes(date,seasonchron,colour = variable))+
  geom_point(shape = 22,aes(fill = variable),colour = "grey20", size = 2)+
  geom_jitter(
    data = EPG[species == "oes"],
    height = .25,aes(date,seasonchron,colour = treatment),alpha = .3)+
  geom_jitter(
    data = EPG[species == "oes"],
    aes(Trtdate,seasonchron),colour = "grey30",fill = "tomato",shape = 21,size = 2,alpha = .2,height = .1)+
  scale_x_date(date_breaks = "1 month")+
  theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.minor = element_blank())

## changing spring 2014
load(file = "data/EPG.rda")
load(file = "data/Treatment.rda")


merged.data <- merge_epg_and_treatment_data(EPG,Treatment,
                                            seasonchron.to.change = "2014-2.spring",new.start = "2014-01-15"
                                            # seasonchron.to.change = NULL
)
EPG <- merged.data$epg.df
Treatment <- merged.data$trt.df

# subset only relevant data -----------------------------------------------

FECR.3w.df <- subset_relevant.data(EPG,n.days.before = 3 * 7,n.days.after = 3 * 7)
FECR.3m.df <- subset_relevant.data(EPG,n.days.before = 12 * 7,n.days.after = 12 * 7)

## keeping 3 weeks around treatment
ggplot(
  melt(change_seasonchron_start.end(seasonchron.to.change = "2014-2.spring",new.start = "2014-01-15"),
       # melt(create_seaonchron.lookup(),
       measure.vars = c("start","end"),value.name = "date"),
  aes(date,seasonchron,colour = variable))+
  geom_point(shape = 22,aes(fill = variable),colour = "grey20", size = 2)+
  geom_jitter(
    data = FECR.3w.df[species == "oes"],
    height = .25,aes(date,seasonchron,colour = period),alpha = .3)+
  geom_jitter(
    data = FECR.3w.df[species == "oes"],
    aes(Trtdate,seasonchron),colour = "grey30",fill = "tomato",shape = 21,size = 2,alpha = .2,height = .1)+
  scale_x_date(date_breaks = "1 month")+
  theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.minor = element_blank())

## keeping 3 months around treatment
ggplot(
  melt(change_seasonchron_start.end(seasonchron.to.change = "2014-2.spring",new.start = "2014-01-15"),
       # melt(create_seaonchron.lookup(),
       measure.vars = c("start","end"),value.name = "date"),
  aes(date,seasonchron,colour = variable))+
  geom_point(shape = 22,aes(fill = variable),colour = "grey20", size = 2)+
  geom_jitter(
    data = FECR.3m.df[species == "oes"],
    height = .25,aes(date,seasonchron,colour = period),alpha = .3)+
  geom_jitter(
    data = FECR.3m.df[species == "oes"],
    aes(Trtdate,seasonchron),colour = "grey30",fill = "tomato",shape = 21,size = 2,alpha = .2,height = .1)+
  scale_x_date(date_breaks = "1 month")+
  theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.minor = element_blank())

## changing spring 2014 & spring 2013
load(file = "data/EPG.rda")
load(file = "data/Treatment.rda")

seasonchron.lookup <- change_seasonchron_start.end(seasonchron.to.change = "2013-2.spring",new.end = "2013-08-20")
seasonchron.lookup <- change_seasonchron_start.end(seasonchron.lookup = seasonchron.lookup,
                                                   seasonchron.to.change = "2014-2.spring",new.start = "2014-01-15")
merged.data <- merge_epg_and_treatment_data(EPG,Treatment,seasonchron.lookup = seasonchron.lookup,
                                            seasonchron.to.change = NULL)
EPG <- merged.data$epg.df
Treatment <- merged.data$trt.df

# subset only relevant data -----------------------------------------------

FECR.3w.df <- subset_relevant.data(EPG,n.days.before = 3 * 7,n.days.after = 3 * 7)
FECR.3m.df <- subset_relevant.data(EPG,n.days.before = 12 * 7,n.days.after = 12 * 7)

## keeping 3 months around treatment
ggplot(
  melt(seasonchron.lookup,
       measure.vars = c("start","end"),value.name = "date"),
  aes(date,seasonchron,colour = variable))+
  geom_point(shape = 22,aes(fill = variable),colour = "grey20", size = 2)+
  geom_jitter(
    data = FECR.3m.df[species == "oes"],
    height = .25,aes(date,seasonchron,colour = period),alpha = .3)+
  geom_jitter(
    data = FECR.3m.df[species == "oes"],
    aes(Trtdate,seasonchron),colour = "grey30",fill = "tomato",shape = 21,size = 2,alpha = .2,height = .1)+
  scale_x_date(date_breaks = "1 month")+
  theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.minor = element_blank())


# plotting averaged epg ---------------------------------------------------

FECR.3w.df.ind <- keep_with_before.after(FECR.3w.df)
FECR.3w.df.ind <- keep_lines_useful_for_FECR.method(FECR.3w.df.ind,"MacIntosh")
FECR.3w.df.ind <- average_per_ind(FECR.3w.df.ind,avg_fun = median)

FECR.3w.df.ind[,.N,by = .(seasonchron,tvnt,period,id,species)]$N
ggplot(FECR.3w.df.ind[species %in% c("oes","trich")],aes(interaction(period,tvnt),epg,colour = period,fill = period))+
  facet_wrap(.~seasonchron + species,nrow = 6,scales = "free")+
  geom_line(aes(group = id),colour = "grey70",alpha = 0.6)+
  geom_point(alpha = 0.3)+
  geom_boxplot(alpha = .4)+
  theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank())
  
# calculate FECR across methods, seasons and parasite species -------------
# only "oes" and "trich" seems to work for some FECR

FECRs.mean.3m <- FECR_from_df.across.season.method(fecr.df = FECR.3m.df[species %in% c("oes","trich")],
                                                   avg_fun = mean,only_with_before.after = TRUE,
                                                   compute.CI = TRUE)
FECRs.mean.3m

FECRs.median.3m <- FECR_from_df.across.season.method(fecr.df = FECR.3m.df[species %in% c("oes","trich")],
                                                     avg_fun = median,only_with_before.after = TRUE,
                                                     compute.CI = TRUE)
FECRs.median.3m

FECRs.mean.3w <- FECR_from_df.across.season.method(fecr.df = FECR.3w.df[species %in% c("oes","trich")],
                                                   avg_fun = mean,only_with_before.after = TRUE,
                                                   compute.CI = TRUE)
FECRs.mean.3w

FECRs.median.3w <- FECR_from_df.across.season.method(fecr.df = FECR.3w.df[species %in% c("oes","trich")],
                                                     avg_fun = median,only_with_before.after = TRUE,
                                                     compute.CI = TRUE)
FECRs.median.3w

# compare results for averaging by mean or median epg (per group or per ind)
FECR.comp <- rbind(
  data.table(FECRs.mean.3m,avg_fun = "mean",n.days = "3m"),
  data.table(FECRs.median.3m,avg_fun = "median",n.days = "3m"),
  data.table(FECRs.mean.3w,avg_fun = "mean",n.days = "3w"),
  data.table(FECRs.median.3w,avg_fun = "median",n.days = "3w")
)

FECR.comp$method <- factor(FECR.comp$method,levels = c("Kochapakdee","Dash","Coles","Cabaret1","Cabaret2","MacIntosh"))

# FECR plots ----------------------------------------------------------
ggplot(FECR.comp,aes(method,FECR,fill = n.days))+
  facet_wrap(seasonchron + avg_fun ~ species,ncol = 4,scales = "free_y")+
  geom_hline(yintercept = 0)+geom_hline(yintercept = 100,lty="dashed",colour="grey50")+
  geom_col(position = position_dodge(.9))+
  geom_errorbar(aes(ymin = FECR.low,ymax = FECR.up),width = 0.2,colour="grey50",position = position_dodge(.9))+
  theme_bw()

ggplot(FECR.comp[avg_fun == "median"],aes(method,FECR,fill = n.days))+
  facet_wrap(seasonchron + avg_fun ~ species,ncol = 2,scales = "free_y")+
  geom_hline(yintercept = 0)+geom_hline(yintercept = 100,lty="dashed",colour="grey50")+
  geom_col(position = position_dodge(.9))+
  geom_errorbar(aes(ymin = FECR.low,ymax = FECR.up),width = 0.2,colour="grey50",position = position_dodge(.9))+
  theme_bw()


