#' Update gene segment status.
#' @description The next state of gene segment status is updated, based on the most current state of other parameters.
#'
#' @return
#' @export
#'
#' @examples
update_gstatus <- function() {

  # parameters
  sum.H <- sum(hh.current)     # sum H at current state
  numer1 <- bb.current * (p1^sum.H) * (1 - p1)^(snp.nums - sum.H)
  numer0 <- (1 - bb.current) * (p0^sum.H) * (1 - p0)^(snp.nums - sum.H)

  # probability
  prob_G <- numer1 / ( numer1 + numer0 )
  newG <- rbinom(1, 1, prob = prob_G)
  ggstatus[iter+1] <<- newG
  assign('gg.current', ggstatus[iter+1], envir = .GlobalEnv)

  return()
}
