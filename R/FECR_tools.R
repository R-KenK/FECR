
# core functions to manipulate data toward FECR wrappers -----------------------

#' Subset lines for which there are data before and after treatment
#' Internal use
#'
#' @param df a data frame to subset
#' @param col.to.keep list of column names for which data presence should be checked for before and after
#' @param col.to.check name of the column stating either "before" or "after". Default is "period".
#' @param bef.aft.names optional if the variable names of the column col.to.check aren't "before" and "after"
#'
#' @return a subsetted data frame
#' @noRd
keep_with_before.after<- function(df,col.to.keep = c("id","period","tvnt","species","seasonchron"),col.to.check = "period",bef.aft.names = c("before","after")){
  is.dt <- data.table::is.data.table(df)
  if (is.dt) {df <- as.data.frame(df);on.exit(setDT(df),add = TRUE,after = FALSE)}
  # pivot table of number of (id,before,after) per tvnt, species and seasonchron
  # transformed into a df with a "Freq" logical column (TRUE when there is at least one line)
  pivot <- check_if_data_per.period(df,col.to.keep)
  
  # adds a "keep" logical column
  pivot <- determine_lines_to_keep(df,pivot,col.to.check = col.to.check,bef.aft.names = bef.aft.names)
  
  # merge the data frame with the pivot (mainly to append the "keep" column)
  df <- merge(df,pivot,by = col.to.keep)
  # subset to only the lines to keep, before removing the "Freq" and "keep" columns
  df <- df[df$keep,][,!(colnames(df) %in% c("Freq","keep"))]
  if (is.dt) {
    setDT(df)
    df
  } else {
    df
  }
}

average_per_group <- function(fecr.df,avg_fun = median) {
  fecr.df[,.(date = mean(date),epg = avg_fun(epg),sd.epg = sd(epg)),by = .(seasonchron,tvnt,period,species)][order(seasonchron)]
}

average_per_ind <- function(fecr.df,avg_fun = median) {
  fecr.df[,.(date=mean(date),epg=median(epg),sd.epg=sd(epg)),by=.(seasonchron,tvnt,id,period,species)][order(seasonchron)]
}

remove_potential_division.by.zero <- function(fecr.df,method = c("Kochapakdee","Dash","Coles","Cabaret1","Cabaret2","MacIntosh")) {
  switch(method,
         "Kochapakdee" = ,
         "Dash" = ,
         "Coles" = fecr.df,
         "Cabaret1" = {
           T1.is.zero <- fecr.df$tvnt == "treated" & fecr.df$period == "before" & fecr.df$epg == 0
           fecr.df <- fecr.df[!T1.is.zero]
           keep_with_before.after(fecr.df)
         },
         "Cabaret2" = ,
         "MacIntosh" = {
           T1.is.zero <- fecr.df$tvnt == "treated" & fecr.df$period == "before" & fecr.df$epg == 0
           C1.is.zero <- fecr.df$tvnt == "control" & fecr.df$period == "before" & fecr.df$epg == 0
           C2.is.zero <- fecr.df$tvnt == "control" & fecr.df$period == "after" & fecr.df$epg == 0
           fecr.df <- fecr.df[!T1.is.zero & !C1.is.zero & !C2.is.zero]
           keep_with_before.after(fecr.df)
         }
  )
}

format_fecr.df <- function(fecr.df,method) {
  fecr.df <- switch(method,
                    "Kochapakdee" = ,
                    "Dash" = ,
                    "Coles" = fecr.df[,c("seasonchron","tvnt","period","species","epg","sd.epg")],
                    "Cabaret1" = ,
                    "Cabaret2" = ,
                    "MacIntosh" = {
                      fecr.df <- fecr.df[,c("seasonchron","tvnt","period","id","species","epg","sd.epg")]
                      fecr.df$id <- factor(fecr.df$id,levels = unique(fecr.df$id))
                      fecr.df[order(seasonchron,tvnt,id,species)]
                    }
  )
  fecr.df$seasonchron <- factor(fecr.df$seasonchron,levels = unique(fecr.df$seasonchron))
  fecr.df$tvnt <- factor(fecr.df$tvnt,levels = c("control","treated"))
  fecr.df$period <- factor(fecr.df$period,levels = c("before","after"))
  fecr.df$species <- factor(fecr.df$species,levels = unique(fecr.df$species))
  fecr.df
}

keep_lines_useful_for_FECR.method <- function(fecr.df,method = c("Kochapakdee","Dash","Coles","Cabaret1","Cabaret2","MacIntosh")) {
  switch(method,
         "Kochapakdee" = ,
         "Cabaret1" = {
           is.T1 <- is.lines_of_group.period(fecr.df,"T1")
           is.T2 <- is.lines_of_group.period(fecr.df,"T2")
           fecr.df[is.T1 | is.T2]
         },
         "Coles" = {
           is.T2 <- is.lines_of_group.period(fecr.df,"T2")
           is.C2 <- is.lines_of_group.period(fecr.df,"C2")
           fecr.df[is.T2 | is.C2]
         },
         "Dash" = ,
         "Cabaret2" = ,
         "MacIntosh" = fecr.df
  )
}


# wrapper for manipulating the data toward FECR wrapper -------------------

prepare_fecr.df <- function(fecr.df,method = c("Kochapakdee","Dash","Coles","Cabaret1","Cabaret2","MacIntosh"),avg_fun = median,only_with_before.after = TRUE) {
  if (only_with_before.after) {
    if (method == "Coles") {
      warning("lines without before and after not removed because this method only considers period == 'after' ")
    } else {
      fecr.df <- keep_with_before.after(fecr.df)
    }
  }
  fecr.df <- keep_lines_useful_for_FECR.method(fecr.df,method)
  fecr.df <- switch(method,
                    "Kochapakdee" = ,
                    "Dash" = ,
                    "Coles" = average_per_group(fecr.df,avg_fun = avg_fun),
                    "Cabaret1" = ,
                    "Cabaret2" = ,
                    "MacIntosh" = average_per_ind(fecr.df,avg_fun = avg_fun)
  )
  fecr.df <- remove_potential_division.by.zero(fecr.df,method = method)
  format_fecr.df(fecr.df,method)
}

# internal functions --------------------------------------------------------

check_if_data_per.period <- function(df,col.to.keep = c("id","period","tvnt","species","seasonchron")) {
  # pivot table of number of (id,before,after) per tvnt, species and seasonchron
  # transformed into a df with the number of line being recorded in the Freq column
  pivot <- data.frame(table(df[,col.to.keep]))
  # turns the number of line into a logical column (TRUE when there is at least one line)
  pivot$Freq <- pivot$Freq > 0
  pivot[,c("seasonchron","tvnt","period","species","id","Freq")]
}

determine_lines_to_keep <- function(df,pivot,col.to.check = "period",bef.aft.names = c("before","after")) {
  # determine which line to keep (when it the individual has data both before and after)
  ## subset to the column "period" == "before" in the default case
  there.is.data_bef <- pivot[pivot[,col.to.check] == bef.aft.names[1],]$Freq
  ## subset to the column "period" == "after" in the default case
  there.is.data_aft <- pivot[pivot[,col.to.check] == bef.aft.names[2],]$Freq
  
  pivot$keep <- NA
  pivot[pivot[,col.to.check] == bef.aft.names[1],]$keep <- there.is.data_bef & there.is.data_aft
  pivot[pivot[,col.to.check] == bef.aft.names[2],]$keep <- there.is.data_bef & there.is.data_aft
  pivot
}


#' Create a random list with a minimum number of individuals
#' Internal use.
#'
#' @param FECR.df data frame of epg data
#' @param s character, season to consider
#' @param group character, group to consider (e.g. "control", "treated")
#' @param per character, period to consider (e.g. "before", "after")
#' @param sp character, parasite species to consider (e.g. "oes", "trich")
#' 
#' @importFrom data.table data.table
#' @importFrom data.table is.data.table
#'
#' @return an ordered vector of individual names
#' @noRd
rand.ind.list <- function(FECR.df,s,group,per,sp,method){
  if (!data.table::is.data.table(FECR.df)) {FECR.df <- data.table::data.table(FECR.df)}
  
  min.ind <- FECR.df[species == sp][,.(.N),by = .(seasonchron,tvnt,period)]
  min.ind <- switch(method,
                    "Cabaret1" = ,
                    "Cabaret2" = min.ind[,.(N = min(N)),by = .(seasonchron,period)],
                    "MacIntosh" = min.ind[,.(N = min(N)),by = .(seasonchron,tvnt)]
  )
  ind.list <- FECR.df[species == sp][,.(.N),by = .(seasonchron,tvnt,period,id)]
  
  sort(sample(ind.list[seasonchron==s&tvnt==group&period==per]$id,min.ind[seasonchron==s&period==per]$N,replace = FALSE))
}

#' Create a random list with a minimum number of individuals
#' Internal use.
#'
#' @param T1 vector of individual pre-treatment epg of the treated group.
#' @param T2 vector of individual post-treatment epg of the treated group.
#' @param C1 vector of individual pre-treatment epg of the control group.
#' @param C2 vector of individual post-treatment epg of the control group.
#'
#' @importFrom data.table data.table
#' @importFrom data.table is.data.table
#'
#' @return a named list (T1,T2,C1,C2) or epg vectors, where there are as many values in C and T
#' @noRd
rand.ind.list.new <- function(T1,T2,C1,C2) {
  if (length(T1) != length(T2) | length(C1) != length(C2)) {
    stop("Not the same number of individuals in before and after for either C or T")
  }
  min.ind <- min(length(T1),length(C1))
  sample(1:min.ind,min.ind)
}


# helper functions --------------------------------------------------------

subset_lines_of_group.period <- function(fecr.df,group.period = c("T1","T2","C1","C2")) {
  fecr.df[is.lines_of_group.period(fecr.df,group.period)]
}

is.lines_of_group.period <- function(fecr.df,group.period = c("T1","T2","C1","C2")) {
  switch(group.period,
         "T1" = fecr.df$tvnt == "treated" & fecr.df$period == "before",
         "T2" = fecr.df$tvnt == "treated" & fecr.df$period == "after",
         "C1" = fecr.df$tvnt == "control" & fecr.df$period == "before",
         "C2" = fecr.df$tvnt == "control" & fecr.df$period == "after"
  )
}

# Wrappers to calculate FECRs ----------------------------------------------------------
FECR_from_df <- function (fecr.df,
                          method = c("Kochapakdee","Dash","Coles","Cabaret1","Cabaret2","MacIntosh"),
                          compute.CI = FALSE,percentile = c(0.025,0.975),boot = 2000) {
  switch(method,
         "Kochapakdee" = {
           T1 <- subset_lines_of_group.period(fecr.df,"T1")$epg
           T2 <- subset_lines_of_group.period(fecr.df,"T2")$epg
           FECR(T1 = T1,T2 = T2,method = method,compute.CI = compute.CI,percentile = percentile,boot = boot)
         },
         "Cabaret1" = {
           T1 <- subset_lines_of_group.period(fecr.df,"T1")$epg
           T2 <- subset_lines_of_group.period(fecr.df,"T2")$epg
           FECR(T1 = T1,T2 = T2,method = method,compute.CI = compute.CI,percentile = percentile,boot = boot)
         },
         "Coles" = {
           T2 <- subset_lines_of_group.period(fecr.df,"T2")$epg
           C2 <- subset_lines_of_group.period(fecr.df,"C2")$epg
           FECR(T2 = T2,C2 = C2,method = method)
         },
         "Dash" = {
           T1 <- subset_lines_of_group.period(fecr.df,"T1")$epg
           T2 <- subset_lines_of_group.period(fecr.df,"T2")$epg
           C1 <- subset_lines_of_group.period(fecr.df,"C1")$epg
           C2 <- subset_lines_of_group.period(fecr.df,"C2")$epg
           
           FECR(T1 = T1,T2 = T2,C1 = C1,C2 = C2,method = method)
         },
         "Cabaret2" = {
           T1 <- subset_lines_of_group.period(fecr.df,"T1")$epg
           T2 <- subset_lines_of_group.period(fecr.df,"T2")$epg
           C1 <- subset_lines_of_group.period(fecr.df,"C1")$epg
           C2 <- subset_lines_of_group.period(fecr.df,"C2")$epg
           
           T.rand.order <- rand.ind.list.new(T1,T2,C1,C2)
           C.rand.order <- rand.ind.list.new(T1,T2,C1,C2)
           T1 <- T1[T.rand.order]
           T2 <- T2[T.rand.order]
           C1 <- C1[C.rand.order]
           C2 <- C2[C.rand.order]
           
           FECR(T1 = T1,T2 = T2,C1 = C1,C2 = C2,method = method,compute.CI = compute.CI,percentile = percentile,boot = boot)
         },
         "MacIntosh" = {
           T1 <- subset_lines_of_group.period(fecr.df,"T1")$epg
           T2 <- subset_lines_of_group.period(fecr.df,"T2")$epg
           C1 <- subset_lines_of_group.period(fecr.df,"C1")$epg
           C2 <- subset_lines_of_group.period(fecr.df,"C2")$epg
           
           FECR(T1 = T1,T2 = T2,C1 = C1,C2 = C2,method = method,compute.CI = compute.CI,percentile = percentile,boot = boot)
         }
  )
}

FECR_from_df.across.season <- function(fecr.df,seasons = unique(fecr.df$seasonchron),
                                       method = c("Kochapakdee","Dash","Coles","Cabaret1","Cabaret2","MacIntosh"),
                                       compute.CI = FALSE,percentile = c(0.025,0.975),boot = 2000) {
  rbind_lapply(
    seasons,
    function(s) {
      rbind_lapply(
        unique(fecr.df$species),
        function(sp) {
          fecr.df.season.sp <- subset(fecr.df,subset = seasonchron == s & species == sp)
          FECR <- FECR_from_df(fecr.df = fecr.df.season.sp,method = method,compute.CI = compute.CI,percentile = percentile,boot = boot)
          if (compute.CI) {
            data.frame(
              seasonchron = s,
              species = sp,
              method = method,
              FECR = FECR,
              FECR.low = ifelse(!is.null(attr(FECR,"CI")[1]),attr(FECR,"CI")[1],NA),
              FECR.up = ifelse(!is.null(attr(FECR,"CI")[2]),attr(FECR,"CI")[2],NA)
            )
          } else {
            data.frame(
              seasonchron = s,
              species = sp,
              method = method,
              FECR = FECR
            )
          }
        }
      )
    }
  )
}

FECR_from_df.across.season.method <- function(fecr.df,avg_fun = mean,only_with_before.after = TRUE,
                                              methods = c("Kochapakdee","Dash","Coles","Cabaret1","Cabaret2","MacIntosh"),
                                              compute.CI = FALSE,percentile = c(0.025,0.975),boot = 2000) {
  rbind_lapply(
    methods,
    function(method) {
      fecr.df <- prepare_fecr.df(fecr.df = fecr.df,method = method,avg_fun = avg_fun,only_with_before.after = only_with_before.after)
      FECR_from_df.across.season(fecr.df = fecr.df,method = method,compute.CI = compute.CI,percentile = percentile,boot = boot)
    }
  )
}
