#' Apply thinning to Markov Chain Monte Carlo simulations.
#'
#' @param sris
#' @param thinby
#'
#' @return
#' @export
#'
#' @examples
thining <- function(sris,thinby=1) {

  if (is.vector(sris)) {

    sris.len <- length(sris)
    tmp <- seq(1,sris.len,by=thinby)
    sris.thin <- sris[tmp]

  } else  if (is.matrix(sris)) {

    sris.len <- dim(sris)[1]
    tmp <- seq(1,sris.len,by=thinby)
    sris.thin <- as.matrix(sris[tmp,])
  }

  return(list(sris.thin))
}

