#' Check if the inputted data are relevant and valid to run the requested FECR method
#'
#' @param T1 pre-treatment mean epg (or vector of individual epg) of the treated group.
#' @param T2 post-treatment mean epg (or vector of individual epg) of the treated group.
#' @param C1 pre-treatment mean epg (or vector of individual epg) of the control group.
#' @param C2 post-treatment mean epg (or vector of individual epg) of the control group.
#' @param method method to base the FECR calculation on:
#' \itemize{
#'  \item{"Kochapakdee"}{the ratio of pre- and post-treatment epg arithmetic means of the treated group}
#'  \item{"Dash"}{the ratio of pre- and post-treatment epg arithmetic means of the treated group, correcting for the evolution of mean epg in the control group.}
#'  \item{"Coles"}{the ratio of the treated and control post-treatment epg arithmetic means}
#'  \item{"Cabaret1"}{the mean ratio of individual pre- and post-treatment epg of the treated group}
#'  \item{"Cabaret2"}{the mean ratio of pre- and post-treatment epg of the treated group, correcting for the individual epg in the control group.}
#'  \item{"MacIntosh"}{Similar to Cabaret2, but update: now run for all pairs of Ts and Cs, to not require Ts and Cs to have the same length.}
#' }
#'
#' @details for more information, cf. Cabaret & Berrag (2004) and the literature they refer to.
#'
#' @return nothing, but should stop or produce a warning when relevant
#' @noRd
check_if_valid_parameters <- function(arg.list) {
  T1 <- arg.list[["T1"]]
  T2 <- arg.list[["T2"]]
  C1 <- arg.list[["C1"]]
  C2 <- arg.list[["C2"]]
  method <- arg.list[["method"]]
  
  switch(method,
         "Kochapakdee" = {
           if (!is.null(C1)|!is.null(C2)) {
             stop("There shouldn't be a control group with this method.")
           }
           if (T1 == 0) {
             stop("T1 cannot be zero.")
           }
         },
         "Dash" = {
           if (any(T1 == 0,C1 == 0,C2 == 0)) {
             stop("T1,C1 or C2 cannot be zero.")
           }
         },
         "Coles" = {
           if (!is.null(T1) | !is.null(C1)) {
             stop("There shouldn't be pretreatment values with this method.")
           }
           if (C2 == 0) {
             stop("C2 cannot be zero.")
           }
         },
         "Cabaret1" = {
           if (!is.null(C1) | !is.null(C2)) {
             stop("There shouldn't be a control group with this method.")
           }
           if (any(length(T1) == 1,length(T2) == 1)) {
             warning("a single value of epg has been provided in place of an expected vector of individual values.")
           }
           if (length(T1) != length(T2)) {
             stop("T1 and T2 are not of the same length.")
           }
         },
         "Cabaret2" = {
           if (any(length(T1) == 1,length(T2) == 1,
                   length(C1) == 1,length(C2) == 1)) {
             warning("a single value of epg has been provided in place of an expected vector of individual values.")
           }
           if (any(length(T1) != length(T2),
                   length(C1) != length(C2),
                   length(T1) != length(C1))) {
             stop("Some of T1,T2,C1,C2 are not of the same length.")
           }
           if (any(T1 == 0,C1 == 0,C2 == 0)) {
             stop("T1,C1 or C2 cannot be zero.")
           }
         },
         "MacIntosh" = {
           if (any(length(T1) == 1,length(T2) == 1,
                   length(C1) == 1,length(C2) == 1)) {
             warning("a single value of epg has been provided in place of an expected vector of individual values.")
           }
           if (any(length(T1) != length(T2),
                   length(C1) != length(C2))) {
             stop("Either T1 vs T2 or C1 vs C2 is not of the same length.")
           }
           if (any(T1 == 0,C1 == 0,C2 == 0)) {
             stop("T1,C1 or C2 cannot be zero.")
           }
         }
  )
}

