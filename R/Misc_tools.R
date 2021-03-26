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

#' Quick optimized equivalent to sample(x,size,replace=TRUE)
#'
#' @param x a vector
#' @param size number of elements to sample
#'
#' @importFrom stats runif
#'
#' @return sample of size "size" taken from x
#' @export
#'
#' @examples
#' quick_sample(1:20,5)
#' # microbenchmark::microbenchmark(
#' #   runif = {(1:20)[ceiling(runif(5,0,20))]},
#' #   quick_sample = quick_sample(1:20,5),
#' #   sample = sample(1:20,5,replace = TRUE),
#' #   times = 1000,control = list("warmup"=100)
#' # )
quick_sample <- function(x,size) {
  x[ceiling(stats::runif(size,0,length(x)))]
}