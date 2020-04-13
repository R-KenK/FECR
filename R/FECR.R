#' Fecal egg count reduction
#' Calculate the FECR according to several methods described (cf. Cabaret & Berrag (2004))
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
#' @details for more informations, cf. Cabaret & Berrag (2004) and the literature they refer to.
#' 
#' @return the fecal egg count reduction, in percentage
#' @export
#'
#' @examples
#' set.seed(42)
#' T1<- rpois(50,100);T1.mean<- mean(T1)
#' T2<- rpois(50,10);T2.mean<- mean(T2)
#' C1<- rpois(50,100);C1.mean<- mean(C1)
#' C2<- rpois(50,90);C2.mean<- mean(C2)
#' 
#' FECR(T1.mean,T2.mean,method = "Kochapakdee")
#' FECR(T1,T2,C1,C2,method = "Dash")
#' FECR(T2 = T2.mean,C2 = C2.mean,method = "Coles")
#' FECR(T2 = T2,C2 = C2,method = "Coles")
#' FECR(T1,T2,method = "Cabaret1")
#' FECR(T1,T2,C1,C2,method = "Cabaret2")
#' FECR(T1,T2,C1,C2,method = "MacIntosh")
#' 
#' T1<- rpois(60,100);T1.mean<- mean(T1)
#' T2<- rpois(60,10);T2.mean<- mean(T2)
#' C1<- rpois(50,100);C1.mean<- mean(C1)
#' C2<- rpois(50,90);C2.mean<- mean(C2)
#' 
#' FECR(T1,T2,C1,C2,method = "MacIntosh")

FECR<- function(T1=NULL,T2=NULL,C1=NULL,C2=NULL,method=c("Kochapakdee","Dash","Coles","Cabaret1","Cabaret2","MacIntosh")){
  method<- match.arg(method)
  switch(method,
         "Kochapakdee" = {
           if(!is.null(C1)|!is.null(C2)){
             stop("There shouldn't be a control group with this method.")
           }
           if(any(length(T1)>1,length(T2)>1)){
             warning("a vector of epg has been provided in place of a mean. The arithmetic mean of the vector has been used.")
             if(length(T1)>1) {T1<- mean(T1)}
             if(length(T2)>1) {T2<- mean(T2)}
           }
           100*(1-T2/T1)
         },
         "Dash" = {
           if(any(length(T1)>1,length(T2)>1,
                  length(C1)>1,length(C2)>1)){
             warning("a vector of epg has been provided in place of a mean. The arithmetic mean of the vector has been used.")
             if(length(T1)>1) {T1<- mean(T1)}
             if(length(T2)>1) {T2<- mean(T2)}
             if(length(C1)>1) {C1<- mean(C1)}
             if(length(C2)>1) {C2<- mean(C2)}
           }
           100*(1-(T2/T1)*(C1/C2))
         },
         "Coles" = {
           if(!is.null(T1)|!is.null(C1)){
             stop("There shouldn't be pretreatment values with this method.")
           }
           if(any(length(T2)>1,length(C2)>1)){
             warning("a vector of epg has been provided in place of a mean. The arithmetic mean of the vector has been used.")
             if(length(T2)>1) {T2<- mean(T2)}
             if(length(C2)>1) {C2<- mean(C2)}
           }
           100*(1-(T2/C2))
         },
         "Cabaret1" = {
           if(!is.null(C1)|!is.null(C2)){
             stop("There shouldn't be a control group with this method.")
           }
           if(any(length(T1)==1,length(T2)==1)){
             warning("a single value of epg has been provided in place of an expected vector of individual values.")
           }
           if(length(T1)!=length(T2)){
             stop("T1 and T2 are not of the same length.")
           }
           n<- length(T1)
           1/n*sum(100*(1-T2/T1))
         },
         "Cabaret2" = {
           if(any(length(T1)==1,length(T2)==1,
                  length(C1)==1,length(C2)==1)){
             warning("a single value of epg has been provided in place of an expected vector of individual values.")
           }
           if(any(length(T1)!=length(T2),
                  length(C1)!=length(C2),
                  length(T1)!=length(C1))){
             stop("Some of T1,T2,C1,C2 are not of the same length.")
           }
           n<- length(T1)
           1/n*sum(100*(1-(T2/T1)*(C1/C2)))
         },
         "MacIntosh" = {
           if(any(length(T1)==1,length(T2)==1,
                  length(C1)==1,length(C2)==1)){
             warning("a single value of epg has been provided in place of an expected vector of individual values.")
           }
           if(any(length(T1)!=length(T2),
                  length(C1)!=length(C2))){
             stop("Either T1 vs T2 or C1 vs C2 is not of the same length.")
           }
           all.pairs<- expand.grid(control=seq_along(C1),treated=seq_along(T1))
           T1<- T1[all.pairs$treated];T2<- T2[all.pairs$treated]
           C1<- C1[all.pairs$control];C2<- C2[all.pairs$control]
           1/nrow(all.pairs)*sum(100*(1-(T2/T1)*(C1/C2)))
         }
  )
}