#' Check if all the individuals of the treated and control groups are present
#'
#' @param T1 a data.frame containing individual epg values of the treated group before treatment,
#'   with one or more value per individual. The data.frame must be in long format (one row per epg
#'   value) and include columns:
#' \itemize{
#'  \item{"id"}{character or factor indicating which individual the epg value is from.}
#'  \item{"epg"}{numeric, the epg value.}
#' }
#' @param T2 a data.frame containing individual epg values of the treated group after treatment,
#'   with one or more value per individual. The data.frame must be in long format (one row per epg
#'   value) and include columns:
#' \itemize{
#'  \item{"id"}{character or factor indicating which individual the epg value is from.}
#'  \item{"epg"}{numeric, the epg value.}
#' }
#' @param C1 a data.frame containing individual epg values of the control group before treatment,
#'   with one or more value per individual. The data.frame must be in long format (one row per epg
#'   value) and include columns:
#' \itemize{
#'  \item{"id"}{character or factor indicating which individual the epg value is from.}
#'  \item{"epg"}{numeric, the epg value.}
#' }
#' @param C2 a data.frame containing individual epg values of the control group after treatment,
#'   with one or more value per individual. The data.frame must be in long format (one row per epg
#'   value) and include columns:
#' \itemize{
#'  \item{"id"}{character or factor indicating which individual the epg value is from.}
#'  \item{"epg"}{numeric, the epg value.}
#' }
#'
#' @return logical, TRUE if all the individuals of the two groups are found in period 1 and 2
#'   (before and after treatment)
#' @noRd
check_ind.presence <- function(T1,T2,C1,C2) {
  setequal(unique(T1$id),unique(T2$id)) & 
    setequal(unique(C1$id),unique(C2$id))
}

#' Gather the group-wise data.frames into one data.table
#'
#' @param T1 a data.frame containing individual epg values of the treated group before treatment,
#'   with one or more value per individual. The data.frame must be in long format (one row per epg
#'   value) and include columns:
#' \itemize{
#'  \item{"id"}{character or factor indicating which individual the epg value is from.}
#'  \item{"epg"}{numeric, the epg value.}
#' }
#' @param T2 a data.frame containing individual epg values of the treated group after treatment,
#'   with one or more value per individual. The data.frame must be in long format (one row per epg
#'   value) and include columns:
#' \itemize{
#'  \item{"id"}{character or factor indicating which individual the epg value is from.}
#'  \item{"epg"}{numeric, the epg value.}
#' }
#' @param C1 a data.frame containing individual epg values of the control group before treatment,
#'   with one or more value per individual. The data.frame must be in long format (one row per epg
#'   value) and include columns:
#' \itemize{
#'  \item{"id"}{character or factor indicating which individual the epg value is from.}
#'  \item{"epg"}{numeric, the epg value.}
#' }
#' @param C2 a data.frame containing individual epg values of the control group after treatment,
#'   with one or more value per individual. The data.frame must be in long format (one row per epg
#'   value) and include columns:
#' \itemize{
#'  \item{"id"}{character or factor indicating which individual the epg value is from.}
#'  \item{"epg"}{numeric, the epg value.}
#' }
#' 
#' @importFrom data.table data.table
#'
#' @return a data.table with columns gp (group.period, one of T1,T2,C1,C2), id, and epg
#' @noRd
rbind_epgs <- function(T1,T2,C1,C2) {
  rbind(
    data.table::data.table(gp = "T1",as.data.frame(T1)[,c("id","epg")]),
    data.table::data.table(gp = "T2",as.data.frame(T2)[,c("id","epg")]),
    data.table::data.table(gp = "C1",as.data.frame(C1)[,c("id","epg")]),
    data.table::data.table(gp = "C2",as.data.frame(C2)[,c("id","epg")])
  )
}

#' Sample randomly only one epg value per individual
#'
#' @param epg.dt a data.table with columns gp (group.period, one of T1,T2,C1,C2), id, and epg,
#'   returned by `rbind_epgs`
#'
#' @return a data.table with columns gp (group.period, one of T1,T2,C1,C2), id, and a single epg
#'   value per id and gp
#' @noRd
sample_per_id <- function(epg.dt) {
  . <- NULL;.SD <- NULL;gp <- NULL;id <- NULL
  epg.dt[,.(epg = .SD$epg %>% sample(1)),by = .(gp,id)]
}

#' Apply method MacIntosh1 to the sampled epg data set
#' Returns a single FECR value, no confidence interval..
#'
#' @param epg.dt a data.table with columns gp (group.period, one of T1,T2,C1,C2), id, and a single
#'   epg value per id and gp
#'
#' @return the FECR value obtained with method MacIntosh1 for the sampled epg data set
#' @noRd
FECR_MacIntosh_Keuk1_from_sample <- function(epg.dt) {
  gp <- NULL
  FECR(epg.dt[gp == 'T1']$epg,
       epg.dt[gp == 'T2']$epg,
       epg.dt[gp == 'C1']$epg,
       epg.dt[gp == 'C2']$epg,
       method = "MacIntosh1",
       compute.CI = FALSE
  )
}
