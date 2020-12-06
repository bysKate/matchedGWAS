#' Update relative risk parameter.
#' @description For J SNPs on the segment, iteratively update relative risk parameter (R_1, ..., R_J). The next state is updated based on the most current state of other parameters.
#'
#' @return
#' @export
#'
#' @examples
update_risk <- function() {

  # by SNP
  for (idx in 1:snp.nums) {

    # Initialize rr.mh.chain with length = (rr.mh.updates+1)
    rr.mh.chain <- c(ri.current[idx], rep(NA,rr.mh.updates))

    # parameters
    hh <- HH[iter+1, idx]

    for (mh.idx in 1:rr.mh.updates ) {

      ri <-  rr.mh.chain[mh.idx]
      # propose
      shape.curr <- 1 + shape.coef * ri
      new.ri <- rgamma(1,shape = shape.curr, rate = rate.coef)

      # rate 1
      tmp1 <-  multinom_ratio(gam = new.ri,
                              alpha = af.current[idx],
                              miu = mr.current[idx],
                              gmodel = gmdl)/
        multinom_ratio(gam = ri,
                       alpha = af.current[idx],
                       miu = mr.current[idx],
                       gmodel = gmdl)
      sumlog <- sum(ssample[idx,]* log(tmp1))
      rate1 <- exp(sumlog)   # product of 9 cells

      # rate 2 - prior
      rate2 <- ifelse(hh,
                      dgamma(new.ri, shape = RR_shape_1[idx], rate = RR_rate_1[idx])/
                        dgamma(ri, shape = RR_shape_1[idx], rate = RR_rate_1[idx])
                      ,
                      dgamma(new.ri, shape = RR_shape_0[idx], rate = RR_rate_0[idx])/
                        dgamma(ri, shape = RR_shape_0[idx], rate = RR_rate_0[idx]) )

      # rate 3 - proposal jump probability
      rate3 <- dgamma(ri,shape = 1 + shape.coef * new.ri, rate = rate.coef) /
        dgamma(new.ri,shape = shape.curr, rate = rate.coef)

      # accept ratio
      rate <- rate1 * rate2 * rate3
      accept <- rbinom(1, 1, min(1, rate))

      # update
      rr.mh.chain[mh.idx+1] <- new.ri * accept + ri * (1-accept)

      # rr.mh.accep_rate <- c(rr.mh.accep_rate, accept)
    }

    # update
    ri.current[idx] <<- rr.mh.chain[rr.mh.updates+1]
    rrisk[iter+1, idx] <<- ri.current[idx]
  }

  return()
}
