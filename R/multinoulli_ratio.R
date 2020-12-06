#' Calculate Multinoulli probabilities.
#'
#' @param miu
#' @param alpha
#' @param gam
#' @param gmodel
#'
#' @return
#' @export
#'
#' @examples
multinoulli_ratio <- function(miu, alpha, gam, gmodel='add0') {
  # gamma distribution(shape, rate)
  # ssample: snp.nums * 9 matrix

  if (gmodel=='add0') {

    increm=(gam-1)/2
    q00 <- (1-alpha)*(1-miu)*(1-alpha)*(1-miu)
    q01 <- (1-alpha)*(1-alpha)*2*miu*(1-miu)* (1+increm) # if R1=gam, then R2 could be negative!
    q02 <- (1-alpha)*(1-alpha)*miu*miu* gam

    q10 <- 2*(1-alpha)*(1-miu)*alpha*miu
    q11 <- 2*alpha*(1-alpha)*(2*miu*miu-2*miu+1)* (1+increm)
    q12 <- 2*alpha*(1-alpha)*miu*(1-miu)* gam

    q20 <- alpha*alpha*miu*miu
    q21 <- alpha*alpha*2*miu*(1-miu)* (1+increm)
    q22 <- alpha*alpha*(1-miu)*(1-miu)* gam

  } else if (gmodel=='rec') {

    q00 <- (1-alpha)*(1-miu)*(1-alpha)*(1-miu)
    q01 <- (1-alpha)*(1-alpha)*2*miu*(1-miu)*1
    q02 <- (1-alpha)*(1-alpha)*miu*miu*gam

    q10 <- 2*(1-alpha)*(1-miu)*alpha*miu
    q11 <- 2*alpha*(1-alpha)*(2*miu*miu-2*miu+1)*1
    q12 <- 2*alpha*(1-alpha)*miu*(1-miu)*gam

    q20 <- alpha*alpha*miu*miu
    q21 <- alpha*alpha*2*miu*(1-miu)* 1
    q22 <- alpha*alpha*(1-miu)*(1-miu)*gam

  } else if (gmodel=='dom') {

    q00 <- (1-alpha)*(1-miu)*(1-alpha)*(1-miu)
    q01 <- (1-alpha)*(1-alpha)*2*miu*(1-miu)*gam
    q02 <- (1-alpha)*(1-alpha)*miu*miu*gam

    q10 <- 2*(1-alpha)*(1-miu)*alpha*miu
    q11 <- 2*alpha*(1-alpha)*(2*miu*miu-2*miu+1)*gam
    q12 <- 2*alpha*(1-alpha)*miu*(1-miu)*gam

    q20 <- alpha*alpha*miu*miu
    q21 <- alpha*alpha*2*miu*(1-miu)*gam
    q22 <- alpha*alpha*(1-miu)*(1-miu)*gam

  } else if (gmodel=='mul') {

    q00 <- (1-alpha)*(1-miu)*(1-alpha)*(1-miu)
    q01 <- (1-alpha)*(1-alpha)*2*miu*(1-miu)* gam
    q02 <- (1-alpha)*(1-alpha)*miu*miu* gam^2

    q10 <- 2*(1-alpha)*(1-miu)*alpha*miu
    q11 <- 2*alpha*(1-alpha)*(2*miu*miu-2*miu+1)* gam
    q12 <- 2*alpha*(1-alpha)*miu*(1-miu)* gam^2

    q20 <- alpha*alpha*miu*miu
    q21 <- alpha*alpha*2*miu*(1-miu)* gam
    q22 <- alpha*alpha*(1-miu)*(1-miu)* gam^2

  } else if (gmodel=='dou') {
    q00 <- (1-alpha)*(1-miu)*(1-alpha)*(1-miu)
    q01 <- (1-alpha)*(1-alpha)*2*miu*(1-miu)*gam
    q02 <- (1-alpha)*(1-alpha)*miu*miu*gam *2

    q10 <- 2*(1-alpha)*(1-miu)*alpha*miu
    q11 <- 2*alpha*(1-alpha)*(2*miu*miu-2*miu+1)*gam
    q12 <- 2*alpha*(1-alpha)*miu*(1-miu)*gam *2

    q20 <- alpha*alpha*miu*miu
    q21 <- alpha*alpha*2*miu*(1-miu)*gam
    q22 <- alpha*alpha*(1-miu)*(1-miu)*gam *2
  }

  deno <- q00+q01+q02+q10+q11+q12+q20+q21+q22
  proportion <- rbind(q00/deno,q01/deno,q02/deno,
                      q10/deno,q11/deno,q12/deno,
                      q20/deno,q21/deno,q22/deno)

  return(proportion)

}
