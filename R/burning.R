#' Apply burning process to Markov Chain Monte Carlo simualtions by disgarding the begining period of sequences.
#'
#' @param sris
#' @param percnt
#'
#' @return
#' @export
#'
#' @examples
burning <- function(sris,percnt) {

  if (is.vector(sris)) {
    sris.len <- length(sris)
    sris.start <- round(sris.len*percnt)

    return(list(sris[sris.start:sris.len]))

  } else {

    sris.len <- dim(sris)[1]
    sris.start <- round(sris.len*percnt)

    return(list(as.matrix(sris[sris.start:sris.len,])))
  }
}
