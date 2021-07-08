#' Average vector of EPG if required
#'
#' @param arg.list list of T1,T2,C1,C2, and method arguments
#'
#' @return list of T1,T2,C1,C2, and method arguments, with vectors of epg averaged if required
#' @noRd
average_if_required <- function(arg.list) {
  T1 <- arg.list[["T1"]]
  T2 <- arg.list[["T2"]]
  C1 <- arg.list[["C1"]]
  C2 <- arg.list[["C2"]]
  method <- arg.list[["method"]]
  
  switch(method,
         "Kochapakdee" = {
           if (any(length(T1) > 1,length(T2) > 1)) {
             warning("a vector of epg has been provided in place of a mean. The arithmetic mean of the vector has been used.")
             if (length(T1) > 1) {T1 <- mean(T1)}
             if (length(T2) > 1) {T2 <- mean(T2)}
           }
           list(T1 = T1,T2 = T2,C1 = C1,C2 = C2,method = method)
         },
         "Dash" = {
           if (
             any(length(T1) > 1,length(T2) > 1,
                 length(C1) > 1,length(C2) > 1
             )
           ) {
             warning("a vector of epg has been provided in place of a mean. The arithmetic mean of the vector has been used.")
             if (length(T1) > 1) {T1 <- mean(T1)}
             if (length(T2) > 1) {T2 <- mean(T2)}
             if (length(C1) > 1) {C1 <- mean(C1)}
             if (length(C2) > 1) {C2 <- mean(C2)}
           }
           list(T1 = T1,T2 = T2,C1 = C1,C2 = C2,method = method)
         },
         "Coles" = {
           if (any(length(T2) > 1,length(C2) > 1)) {
             warning("a vector of epg has been provided in place of a mean. The arithmetic mean of the vector has been used.")
             if (length(T2) > 1) {T2 <- mean(T2)}
             if (length(C2) > 1) {C2 <- mean(C2)}
           }
           list(T1 = T1,T2 = T2,C1 = C1,C2 = C2,method = method)
         },
         "Cabaret1" = ,
         "Cabaret2" = ,
         "MacIntosh" = list(T1 = T1,T2 = T2,C1 = C1,C2 = C2,method = method)
    )
}