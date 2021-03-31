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
#' @details for more information, cf. Cabaret & Berrag (2004) and the literature they refer to.
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

FECR <- function(T1 = NULL,T2 = NULL,C1 = NULL,C2 = NULL,
                 method = c("Kochapakdee","Dash","Coles","Cabaret1","Cabaret2","MacIntosh"),
                 compute.CI = FALSE,percentile = c(0.025,0.975),boot = 2000,
                 boot.original.data = FALSE) {
  method <- match.arg(method)
  check_if_valid_parameters(T1 = T1,T2 = T2,C1 = C1,C2 = C2,method = method)
  
  switch(method,
         "Kochapakdee" = FECR_Kochapakdee(T1,T2),
         "Dash" = FECR_Dash(T1,T2,C1,C2),
         "Coles" = FECR_Coles(T2,C2),
         "Cabaret1" = FECR_Cabaret1(T1,T2,compute.CI = compute.CI,percentile = percentile,boot = boot,boot.original.data = boot.original.data),
         "Cabaret2" = FECR_Cabaret2(T1,T2,C1,C2,compute.CI = compute.CI,percentile = percentile,boot = boot,boot.original.data = boot.original.data),
         "MacIntosh" = FECR_MacIntosh_Keuk(T1,T2,C1,C2,compute.CI = compute.CI,percentile = percentile,boot = boot,boot.original.data = boot.original.data)
  )
}

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
check_if_valid_parameters <- function(T1 = NULL,T2 = NULL,
                                      C1 = NULL,C2 = NULL,
                                      method) {
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

#' Calculate FECR according the Kochapakdee method
#'
#' @param T1 pre-treatment mean epg of the treated group.
#' @param T2 post-treatment mean epg of the treated group.
#'
#' @details for more informations, cf. Cabaret & Berrag (2004) and the literature they refer to.
#'
#' @return the FECR with this method
#' @noRd
FECR_Kochapakdee <- function(T1,T2) {
  if (any(length(T1) > 1,length(T2) > 1)) {
    warning("a vector of epg has been provided in place of a mean. The arithmetic mean of the vector has been used.")
    if (length(T1) > 1) {T1 <- mean(T1)}
    if (length(T2) > 1) {T2 <- mean(T2)}
  }
  r <- 100 * (1 - T2/T1)
  ifelse(r >= 0,r,0) # truncates at zero (when T2 > T1)
}

#' Calculate FECR according the Dash method
#'
#' @param T1 pre-treatment mean epg of the treated group.
#' @param T2 post-treatment mean epg of the treated group.
#'
#' @details for more information, cf. Cabaret & Berrag (2004) and the literature they refer to.
#'
#' @return the FECR with this method
#' @noRd
FECR_Dash <- function(T1,T2,C1,C2) {
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
  r <- 100 * (1 - (T2/T1) / (C2/C1))
  ifelse(r >= 0,r,0)
}

#' Calculate FECR according the Coles method
#'
#' @param T1 pre-treatment mean epg of the treated group.
#' @param T2 post-treatment mean epg of the treated group.
#'
#' @details for more information, cf. Cabaret & Berrag (2004) and the literature they refer to.
#'
#' @return the FECR with this method
#' @noRd
FECR_Coles <- function(T2,C2) {
  if (any(length(T2) > 1,length(C2) > 1)) {
    warning("a vector of epg has been provided in place of a mean. The arithmetic mean of the vector has been used.")
    if (length(T2) > 1) {T2 <- mean(T2)}
    if (length(C2) > 1) {C2 <- mean(C2)}
  }
  r <- 100 * (1 - T2/C2)
  ifelse(r >= 0,r,0)
}


#' function containing the formula for the Cabaret1 method
#'
#' @param n integer, number of individuals in both treated and control group
#' @param T1 vector of individual pre-treatment epg of the treated group.
#' @param T2 vector of individual post-treatment epg of the treated group.
#'
#' @return the FECR with this method
#' @noRd
Cabaret1_fun <- function(n,T1,T2) {
  r <- mean(100 * (1 - T2/T1))
  ifelse(r >= 0,r,0)
}

#' Calculate FECR according the Cabaret1 method
#'
#' @param T1 vector of individual pre-treatment epg of the treated group.
#' @param T2 vector of individual post-treatment epg of the treated group.
#' @param compute.CI logical. Should confidence interval be calculated? 
#' @param percentile quantiles to use if compute.CI is TRUE, unused otherwise.
#' @param boot number of bootstrap iteration to compute if compute.CI is TRUE, unused otherwise.
#'
#' @details for more information, cf. Cabaret & Berrag (2004) and the literature they refer to.
#'
#' @return
#' @noRd
FECR_Cabaret1 <- function(T1,T2,
                          compute.CI = FALSE,
                          percentile = c(0.025,0.975),
                          boot = 2000,
                          boot.original.data = FALSE) {
  n <- length(T1)
  FECR <- Cabaret1_fun(n,T1,T2)
  if (compute.CI) {
    if (boot.original.data) {
      bootstrap <- pbapply::pbsapply(
        1:boot,
        function(b) {
          T1.resampled <- quick_sample(T1,n)
          T2.resampled <- quick_sample(T1,n)
          Cabaret1_fun(n,T1.resampled,T2.resampled)
        }
      )
    } else {
      T2.T1 <- T2 / T1
      bootstrap <- pbapply::pbsapply(
        1:boot,
        function(b) {
          T2.T1.resampled <- quick_sample(T2.T1,n)
          r <- mean(100 * (1 - T2.T1.resampled))
          ifelse(r >= 0,r,0)
        }
      )
    }
    attr(FECR,"CI") <- c(
      lower = stats::quantile(bootstrap,percentile[1]),
      upper = stats::quantile(bootstrap,percentile[2])
    )
  }
  FECR
}

#' function containing the formula for the Cabaret2 method
#'
#' @param n integer, number of individuals in both treated and control group
#' @param T1 vector of individual pre-treatment epg of the treated group.
#' @param T2 vector of individual post-treatment epg of the treated group.
#' @param C1 vector of individual pre-treatment epg of the control group.
#' @param C2 vector of individual post-treatment epg of the control group.
#'
#' @return the FECR with this method
#' @noRd
Cabaret2_fun <- function(n,T1,T2,C1,C2) {
  r <- mean(100 * (1 - (T2/T1) / (C2/C1)))
  ifelse(r >= 0,r,0)
}

#' Calculate FECR according the Cabaret2 method
#'
#' @param T1 vector of individual pre-treatment epg of the treated group.
#' @param T2 vector of individual post-treatment epg of the treated group.
#' @param C1 vector of individual pre-treatment epg of the control group.
#' @param C2 vector of individual post-treatment epg of the control group.
#' @param compute.CI logical. Should confidence interval be calculated? 
#' @param percentile quantiles to use if compute.CI is TRUE, unused otherwise.
#' @param boot number of bootstrap iteration to compute if compute.CI is TRUE, unused otherwise.
#'
#' @details for more information, cf. Cabaret & Berrag (2004) and the literature they refer to.
#'
#' @return the FECR with this method, with or without its bootstrapped confidence interval
#' @noRd
FECR_Cabaret2 <- function(T1,T2,C1,C2,
                          compute.CI = FALSE,
                          percentile = c(0.025,0.975),
                          boot = 2000,
                          boot.original.data = FALSE) {
  n <- length(T1)
  
  # # to produce a random pairing of C and T, but this is not the behavior of
  # # their original equation:
  # C.order <- sample(1:n,n);T.order <- sample(1:n,n)
  # C1 <- C1[C.order];C2 <- C2[C.order]
  # T1 <- T1[T.order];T2 <- T2[T.order]
  
  FECR <- Cabaret2_fun(n,T1,T2,C1,C2)
  
  if (compute.CI) {
    if (boot.original.data) {
      T2.T1 <- T2 / T1
      C2.C1 <- C2 / C1
      bootstrap <- pbapply::pbsapply(
        1:boot,
        function(b) {
          T.resampled <- quick_sample(T2.T1,n)
          C.resampled <- quick_sample(C2.C1,n)
          r <- mean(100 * (1 - T.resampled / C.resampled))
          ifelse(r >= 0,r,0)
        }
      )
    } else {
      TC <- (T2 / T1) / (C2 / C1)
      bootstrap <- pbapply::pbsapply(
        1:boot,
        function(b) {
          TC.resampled <- quick_sample(TC,n)
          r <- mean(100 * (1 - TC.resampled))
          ifelse(r >= 0,r,0)
        }
      )
    }
    attr(FECR,"CI") <- c(
      stats::quantile(bootstrap,percentile[1]),
      stats::quantile(bootstrap,percentile[2])
    )
  }
  FECR
}

#' function containing the formula for the modified Cabaret2 (MacIntosh & Keuk) method
#'
#' @param n.m integer, product of the number of individuals in treated and control group (e.g. for respective group size n and m, n.m = n * m)
#' @param T1 vector of individual pre-treatment epg of the treated group.
#' @param T2 vector of individual post-treatment epg of the treated group.
#' @param C1 vector of individual pre-treatment epg of the control group.
#' @param C2 vector of individual post-treatment epg of the control group.
#'
#' @return the FECR with this method
#' @noRd
MacIntosh_Keuk_fun <- function(T.all.pairs,C.all.pairs) {
  r <- mean(100 * (1 - T.all.pairs / C.all.pairs))
  # r
  ifelse(r >= 0,r,0)
}

#' Calculate FECR according the modified Cabaret2 (MacIntosh & Keuk method)
#'
#' @param T1 vector of individual pre-treatment epg of the treated group.
#' @param T2 vector of individual post-treatment epg of the treated group.
#' @param C1 vector of individual pre-treatment epg of the control group.
#' @param C2 vector of individual post-treatment epg of the control group.
#' @param compute.CI logical. Should confidence interval be calculated? 
#' @param percentile quantiles to use if compute.CI is TRUE, unused otherwise.
#' @param boot number of bootstrap iteration to compute if compute.CI is TRUE, unused otherwise.
#'
#' @details for more information, cf. Cabaret & Berrag (2004) and the literature they refer to.
#'
#' @return the FECR with this method, with or without its bootstrapped confidence interval
#' @noRd
FECR_MacIntosh_Keuk <- function(T1,T2,C1,C2,
                          compute.CI = FALSE,
                          percentile = c(0.025,0.975),
                          boot = 2000,
                          boot.original.data = FALSE) {
  # list all pairings in a data frame
  all.pairs <- expand.grid(control = seq_along(C1),treated = seq_along(T1))
  n.m <- nrow(all.pairs)
  
  # reshape into longer vectors (length = n² (or = n * m = length(C) * length(T) if they don't have the same length)),
  # covering all T and C pairs when aligned
  T.all.pairs <- T2[all.pairs$treated] / T1[all.pairs$treated]
  C.all.pairs <- C2[all.pairs$control] / C1[all.pairs$control]
  
  FECR <- MacIntosh_Keuk_fun(T.all.pairs,C.all.pairs)
  if (compute.CI) {
    if (boot.original.data) {
      T2.T1 <- T2 / T1
      C2.C1 <- C2 / C1
      bootstrap <- pbapply::pbsapply(
        1:boot,
        function(b) {
          T.resampled <- quick_sample(T2.T1,length(T2.T1))
          C.resampled <- quick_sample(C2.C1,length(C2.C1))
          
          T.all.pairs.resampled <- T.resampled[all.pairs$treated]
          C.all.pairs.resampled <- C.resampled[all.pairs$control]
          TC.all.pairs.resampled <- T.all.pairs.resampled / C.all.pairs.resampled
          
          r <- mean(100 * (1 - TC.all.pairs.resampled))
          ifelse(r >= 0,r,0)
        }
      )
    } else {
      TC.all.pairs <- T.all.pairs / C.all.pairs
      bootstrap <- pbapply::pbsapply(
        1:boot,
        function(b){
          TC.all.pairs.resampled <- quick_sample(TC.all.pairs,n.m)
          r <- mean(100 * (1 - TC.all.pairs.resampled))
          ifelse(r >= 0,r,0)
        }
      )
    }
    attr(FECR,"CI") <- c(
      stats::quantile(bootstrap,percentile[1]),
      stats::quantile(bootstrap,percentile[2])
    )
  }
  FECR
}
