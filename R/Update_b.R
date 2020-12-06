#' Update probability parameter b.
#' @description The next state of parameter b is updated, based on the most current state of other parameters.
#'
#' @return
#' @export
#'
#' @examples
update_b <- function() {
  # posterior: s1=G + alpha   s2=beta-G+1
  alpha_new <- ggstatus[iter] + alpha_b
  beta_new <- beta_b - ggstatus[iter] + 1

  bb[iter+1] <<- rbeta(1, shape1 = alpha_new, shape2 = beta_new)

  assign('bb.current', bb[iter+1], envir = .GlobalEnv)
  return()
}
