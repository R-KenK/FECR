#' TO WRITE
#'
#' @param seasonchron TO WRITE
#' @param start TO WRITE
#' @param end TO WRITE
#'
#' @return TO WRITE
#' @noRd
create_seaonchron.lookup <- function(
  seasonchron = c(
    "2012-4.fall","2013-1.winter","2013-2.spring","2013-3.summer",
    "2013-4.fall","2014-1.winter","2014-2.spring","2014-3.summer",
    "2014-4.fall","2015-2.spring","2015-3.summer","2015-4.fall"
  ),
  start = as.Date(
    c(
      "2012-09-01",
      "2013-01-01",
      "2013-04-01",
      "2013-07-01", # some epg might be included to "2013-2.spring" if this one is later PLUs  # this one doesn't match the later years
      "2013-09-01",
      "2014-01-01",
      "2014-04-01", # some epg might be included to "2014-1.winter" if this one is later (i.e. ~ 2014-05-15)
      "2014-06-01", # this one doesn't match the previous year
      "2014-09-01",
      "2015-01-01",
      "2015-06-01",
      "2015-09-01"
    )
  ),
  end = as.Date(
    c(
      "2012-12-31",
      "2013-03-31",
      "2013-06-30",
      "2013-08-31",
      "2013-12-31",
      "2014-03-31",
      "2014-05-31",
      "2014-08-31",
      "2014-12-21",
      "2015-05-31",
      "2015-08-31",
      "2015-12-31"
    )
  )) {
  data.table(
    seasonchron = seasonchron,
    start = start,
    end = end
  )
}

#' TO WRITE
#'
#' @param seasonchron.lookup TO WRITE
#' @param seasonchron.to.change TO WRITE
#' @param new.start TO WRITE
#' @param new.end TO WRITE
#'
#' @return TO WRITE 
#' @noRd
change_seasonchron_start.end <- function(seasonchron.lookup = create_seaonchron.lookup(),seasonchron.to.change = NULL,new.start = NULL,new.end = NULL) {
  # returns the inputted seasonchron.lookup if no changes are requested
  if (is.null(seasonchron.to.change)) {return(seasonchron.lookup)}
  
  to.change.index <- which(seasonchron.lookup$seasonchron == seasonchron.to.change)
  if (is.null(new.start)) {
    new.start <- seasonchron.lookup[to.change.index,]$start
  }
  if (is.null(new.end)) {
    new.end <- seasonchron.lookup[to.change.index,]$end
  }
  new.start <- as.Date(new.start)
  new.end <- as.Date(new.end)
  
  seasonchron.lookup[to.change.index,]$start <- new.start
  seasonchron.lookup[to.change.index,]$end <- new.end
  
  if (to.change.index - 1 > 0) {
    previous.season <- seasonchron.lookup$seasonchron[to.change.index - 1]
    seasonchron.lookup[seasonchron.lookup$seasonchron == previous.season,]$end <- new.start - 1
  }
  if (to.change.index + 1 <= nrow(seasonchron.lookup)) {
    next.season <- seasonchron.lookup$seasonchron[to.change.index + 1]
    seasonchron.lookup[seasonchron.lookup$seasonchron == next.season,]$start <- new.end + 1
  }
  seasonchron.lookup
}

#' TO WRITE
#'
#' @param seasonchron.lookup TO WRITE
#' @param date.vec TO WRITE
#'
#' @return TO WRITE 
#' @noRd
find_seasonchron.index <- function(seasonchron.lookup,date.vec) {
  sapply(
    date.vec,
    function(d) {
      is.in.between <- d >= seasonchron.lookup$start & d <= seasonchron.lookup$end
      if (all(!is.in.between) || table(is.in.between)["TRUE"] != 1) {stop("No season or multiple seasons found for a sample")}
      which(is.in.between)
    }
  )
}

#' TO WRITE
#'
#' @param epg.df TO WRITE
#' @param seasonchron.lookup TO WRITE
#'
#' @return
#' @export
#'
#' @examples
match_epg_seasonchron <- function(epg.df,seasonchron.lookup) {
  original.col.order <- colnames(epg.df)
  # gets a single (averaged) date per sample
  to.match <- epg.df[,.(date = as.Date(mean(date)),sd = sd(date)),by = (sample)]
  if (any(to.match$sd != 0)) {stop("At least one sample has more than one date across parasite species.")}
  
  seasonchron.index <- find_seasonchron.index(seasonchron.lookup,to.match$date)
  matched <- cbind(to.match,seasonchron = seasonchron.lookup$seasonchron[seasonchron.index])
  epg.df <- merge(epg.df[,-c("seasonchron")],matched[,c("sample","seasonchron")],by = "sample")
  epg.df[,original.col.order,with = FALSE]
}


#' TO WRITE
#'
#' @param Treatment.df TO WRITE
#' @param seasonchron.lookup TO WRITE
#'
#' @return TO WRITE 
#' @noRd
match_treatment_seasonchron <- function(Treatment.df,seasonchron.lookup) {
  original.col.order <- colnames(Treatment.df)
  seasonchron.index <- find_seasonchron.index(seasonchron.lookup,Treatment.df$Trtdate)
  Treatment.df <- cbind(subset(Treatment.df,select = setdiff(colnames(Treatment.df),"seasonchron")),seasonchron = seasonchron.lookup$seasonchron[seasonchron.index])
  Treatment.df[,original.col.order,with = FALSE]
}

#' TO WRITE
#'
#' @param epg.df TO WRITE
#' @param trt.df TO WRITE
#'
#' @return TO WRITE 
#' @noRd
keep_with_epg_and_treatment <- function(epg.df,trt.df) {
  id.or.control <- ifelse(epg.df$tvnt == "treated",epg.df$id,"control")
  there.is.date<- sapply(seq_len(nrow(epg.df)),
                         function(i) {
                           seasonchron.epg<- epg.df$seasonchron[i];
                           id.epg<- id.or.control[i];
                           treatment<- epg.df$treatment[i];D<- epg.df$date[i];
                           
                           D.list<- trt.df[trt.df$id == id.epg & trt.df$seasonchron == seasonchron.epg,]$Trtdate
                           cat(paste0("\rcalculating date ",i,"/",nrow(epg.df)))
                           length(D.list)!=0
                         }
  )
  cat("\n")
  epg.df[there.is.date,]
}

#' TO WRITE
#'
#' @param epg.df TO WRITE
#' @param trt.df TO WRITE
#'
#' @return TO WRITE 
#' @noRd
match_epg_with_treatment.date <- function(epg.df,trt.df) {
  epg.df$Trtdate <- as.Date(
    sapply(seq_len(nrow(epg.df)),
           function(i) {
             seasonchron.epg <- epg.df$seasonchron[i];
             id.epg <- epg.df$id[i];tvnt.epg <- epg.df$tvnt[i];
             treatment <- epg.df$treatment[i];D <- epg.df$date[i];
             
             D.list <- max(as.Date(trt.df[trt.df$id == ifelse(tvnt.epg == "treated",id.epg,"control") & trt.df$seasonchron == seasonchron.epg,]$Trtdate))
             cat(paste0("\rcalculating date ",i,"/",nrow(epg.df)))
             closest.date(D,D.list)
           }
    ),
    origin="1970-01-01"
  )
  cat("\n")
  epg.df
}

#' TO WRITE
#'
#' @param epg.df TO WRITE
#' @param trt.df TO WRITE
#' @param seasonchron.lookup TO WRITE
#' @param seasonchron.to.change TO WRITE
#' @param new.start TO WRITE
#' @param new.end TO WRITE
#'
#' @return TO WRITE 
#' @noRd
merge_epg_and_treatment_data <- function(epg.df,trt.df,
                                         seasonchron.lookup = create_seaonchron.lookup(),
                                         seasonchron.to.change = NULL,
                                         new.start = NULL,
                                         new.end = NULL) {
  setDT(epg.df)
  setDT(trt.df)
  
  # change seasonchron date as requested (default is nothing changed)
  seasonchron.lookup <- change_seasonchron_start.end(seasonchron.lookup,seasonchron.to.change,new.start,new.end)
  
  
  # match samples to the season(chron) defined in seasonchron.lookup
  epg.df <- match_epg_seasonchron(epg.df,seasonchron.lookup)
  
  # match samples to the season(chron) defined in seasonchron.lookup
  trt.df <- match_treatment_seasonchron(trt.df,seasonchron.lookup)
  
  # Identify and subset epg data with treatment for the specific season
  epg.df <- keep_with_epg_and_treatment(epg.df,trt.df)
  # Match epg data with the relevant individual or average treatment date
  epg.df <- match_epg_with_treatment.date(epg.df,trt.df)
  list(epg.df = epg.df,trt.df = trt.df)
}
