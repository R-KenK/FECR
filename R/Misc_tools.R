#' Determine the closest date
#'
#' @param D a date
#' @param D.list a vector of dates
#'
#' @return the closest date
#' @export
#'
#' @examples
#' closest.date("1991-03-30",D.list=as.Date(runif(10,0,20),origin="1991-03-15"))
closest.date<- function(D,D.list){
  D.list[which.min(as.Date(D)-as.Date(D.list))]
}

#' Row bind list of data frames
#' wrapper to one-function do.call rbind over a lapply list
#'
#' @param X a list. See details \link[base]{lapply}.
#' @param FUN a function to subset data frames (or data tables). See details \link[base]{lapply}.
#'
#' @return a row bound data frame
#' @export
#'
#' @examples
#' set.seed(42)
#'
#' X<- lapply(1:3,function(i) list(int = 42,df = data.frame(x = runif(10,0,1),y = runif(10,0,1))))
#' rbind_lapply(X,function(x) x$df)
rbind_lapply<- function(X,FUN){
  do.call(rbind,lapply(X = X,FUN = FUN))
}

#' Subset lines for which there are data before and after treatment
#' Internal use
#'
#' @param df a data frame to subset
#' @param col.to.keep list of column names for which data presence should be checked for before and after
#' @param col.to.check name of the colum stating either "before" or "after". Default is "period".
#' @param bef.aft.names optional if the variable names of the colum col.to.check aren't "before" and "after"
#'
#' @return a subsetted data frame
#' @export
#'
#' @examples
#' #internal use
keep_with_before.after<- function(df,col.to.keep=c("id","period","tvnt","species","seasonchron"),col.to.check="period",bef.aft.names=c("before","after")){
  pivot<- data.frame(table(df[,col.to.keep]))
  pivot$Freq<- pivot$Freq>0
  there.is.data_bef<- pivot[pivot[,col.to.check]==bef.aft.names[1],]$Freq
  there.is.data_aft<- pivot[pivot[,col.to.check]==bef.aft.names[2],]$Freq
  
  pivot$keep<-NA
  pivot[pivot[,col.to.check]==bef.aft.names[1],]$keep<- there.is.data_bef&there.is.data_aft
  pivot[pivot[,col.to.check]==bef.aft.names[2],]$keep<- there.is.data_bef&there.is.data_aft
  
  df<- merge(df,pivot,by = col.to.keep)
  df[df$keep,][,!(colnames(df) %in% c("Freq","keep"))]
}

#' Create a random list with a minimum number of individuals
#' Internal use.
#'
#' @param FECR.ind data frame of epg data
#' @param s character, season to consider
#' @param group character, group to consider (e.g. "control", "treated")
#' @param per character, period to consider (e.g. "before", "after")
#' @param sp character, parasite species to consider (e.g. "oes", "trich")
#' 
#' @importFrom data.table data.table
#' @importFrom data.table is.data.table
#'
#' @return an ordered vector of inedividual names
#' @export
#'
#' @examples
#' #internal use.
rand.ind.list<- function(FECR.ind,s,group,per,sp){
  if(!data.table::is.data.table(FECR.ind)){FECR.ind<- data.table::data.table(FECR.ind)}
  min.ind<- FECR.ind[species==sp][,.(.N),by=.(seasonchron,tvnt,period)][,.(N=min(N)),by=.(seasonchron,period)]
  ind.list<- FECR.ind[species==sp][,.(.N),by=.(seasonchron,tvnt,period,id)]
  
  sort(sample(ind.list[seasonchron==s&tvnt==group&period==per]$id,min.ind[seasonchron==s&period==per]$N,replace = FALSE))
}
