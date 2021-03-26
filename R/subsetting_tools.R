#' TO WRITE
#'
#' @param epg.df TO WRITE
#' @param n.days.before TO WRITE
#' @param n.days.after TO WRITE
#'
#' @return TO WRITE
#' @noRd
keep_within_time.window <- function(epg.df,n.days.before = 3 * 7,n.days.after = 3 * 7) {
  EPG[EPG$date >= (EPG$Trtdate - n.days.before) & EPG$date <= (EPG$Trtdate + n.days.after),]   # confirm timing for inclusion
}

#' TO WRITE
#'
#' @param epg.df TO WRITE
#'
#' @return TO WRITE
#' @noRd
determine_before.after <- function(epg.df) {
  epg.df$period <- ifelse(epg.df$date <= epg.df$Trtdate,"before","after")   # what about on the day => "before" to confirm.
  epg.df$period <- factor(epg.df$period,levels = c("before","after"))
  epg.df
}


#' TO WRITE
#'
#' @param epg.df TO WRITE
#' @param n.days.before TO WRITE
#' @param n.days.after TO WRITE
#'
#' @return TO WRITE
#' @noRd
subset_relevant.data <- function(epg.df,n.days.before = 3 * 7,n.days.after = 3 * 7) {
  # Keep only epg data when the collection date is in the relevant time window with respect to treatment date
  epg.df <- keep_within_time.window(epg.df,n.days.before,n.days.after)
  
  # Determine before and after periods based on matched treatment date (revamping the previous sort-of-manual period attribution)
  epg.df <- determine_before.after(epg.df)
  epg.df[,c("date","season","seasonchron","sample","id","tvnt","Trtdate","period","species","epg")]
}
