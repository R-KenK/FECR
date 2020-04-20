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
#' @param compute.CI logical. Should confidence interval be calculated? 
#' @param percentile quantiles to use if compute.CI is TRUE, unused otherwise.
#' @param boot number of bootstrap iteration to compute if compute.CI is TRUE, unused otherwise.
#'
#' @details for more informations, cf. Cabaret & Berrag (2004) and the literature they refer to.
#' 
#' @return the fecal egg count reduction, in percentage (if compute.CI is TRUE, a vector of lower and upper limit quantiles are passed as a "CI" attribute).
#' @export
#' 
#' @importFrom pbapply pbsapply
#' @importFrom stats quantile
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
#' FECR(T1,T2,C1,C2,method = "Cabaret2",compute.CI=TRUE)
#' FECR(T1,T2,C1,C2,method = "MacIntosh",compute.CI=TRUE)
#' 
#' T1<- rpois(60,100);T1.mean<- mean(T1)
#' T2<- rpois(60,10);T2.mean<- mean(T2)
#' C1<- rpois(50,100);C1.mean<- mean(C1)
#' C2<- rpois(50,90);C2.mean<- mean(C2)
#' 
#' FECR(T1,T2,C1,C2,method = "MacIntosh")
#' FECR(T1,T2,C1,C2,method = "MacIntosh",compute.CI=TRUE)

FECR<- function(T1=NULL,T2=NULL,C1=NULL,C2=NULL,method=c("Kochapakdee","Dash","Coles","Cabaret1","Cabaret2","MacIntosh"),compute.CI=FALSE,percentile=c(0.025,0.975),boot=2000){
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
           r<- 100*(1-T2/T1)
           ifelse(r>=0,r,0)
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
           r<- 100*(1-(T2/T1)*(C1/C2))
           ifelse(r>=0,r,0)
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
           r<- 100*(1-(T2/C2))
           ifelse(r>=0,r,0)
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
           r<- 1/n*sum(100*(1-T2/T1))
           if(compute.CI){
             bootstrap<- pbapply::pbsapply(
               1:boot,
               function(b){
                 T1<- sample(T1,n,replace = TRUE);T2<- sample(T2,n,replace = TRUE);
                 r<- 1/n*sum(100*(1-T2/T1))
                 ifelse(r>=0,r,0)
               }
             )
             FECR.mean<- ifelse(r>=0,r,0);attr(FECR.mean,"CI")<-c(lower = stats::quantile(bootstrap,percentile[1]),
                                                                  upper = stats::quantile(bootstrap,percentile[2]))
             FECR.mean
           }else{
             ifelse(r>=0,r,0)
           }
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
           r<- 1/n*sum(100*(1-(T2/T1)*(C1/C2)))
           ifelse(r>=0,r,0)
           if(compute.CI){
             bootstrap<- pbapply::pbsapply(
               1:boot,
               function(b){
                 C1<- sample(C1,n,replace = TRUE);C2<- sample(C2,n,replace = TRUE);
                 T1<- sample(T1,n,replace = TRUE);T2<- sample(T2,n,replace = TRUE);
                 r<- 1/n*sum(100*(1-(T2/T1)*(C1/C2)))
                 ifelse(r>=0,r,0)
               }
             )
             FECR.mean<- ifelse(r>=0,r,0);attr(FECR.mean,"CI")<-c(stats::quantile(bootstrap,percentile[1]),
                                                                  stats::quantile(bootstrap,percentile[2]))
             FECR.mean
           }else{
             ifelse(r>=0,r,0)
           }
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
           r<- 1/nrow(all.pairs)*sum(100*(1-(T2/T1)*(C1/C2)))
           ifelse(r>=0,r,0)
           if(compute.CI){
             bootstrap<- pbapply::pbsapply(
               1:boot,
               function(b){
                 C1<- sample(C1,length(C1),replace = TRUE);C2<- sample(C2,length(C2),replace = TRUE);
                 T1<- sample(T1,length(T1),replace = TRUE);T2<- sample(T2,length(T2),replace = TRUE);
                 
                 T1<- T1[all.pairs$treated];T2<- T2[all.pairs$treated]
                 C1<- C1[all.pairs$control];C2<- C2[all.pairs$control]
                 r<- 1/nrow(all.pairs)*sum(100*(1-(T2/T1)*(C1/C2)))
                 ifelse(r>=0,r,0)
               }
             )
             FECR.mean<- ifelse(r>=0,r,0);attr(FECR.mean,"CI")<-c(stats::quantile(bootstrap,percentile[1]),
                                                                  stats::quantile(bootstrap,percentile[2]))
             FECR.mean
             # bootstrap
           }else{
             ifelse(r>=0,r,0)
           }
         }
  )
}