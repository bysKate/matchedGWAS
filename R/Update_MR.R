#' Update mutation rate parameter.
#'
#' @description For J SNPs on the segment, iteratively update relative risk parameter (M_1, ..., M_J). The next state is updated based on the most current state of other parameters.
#' @return
#' @export
#'
#' @examples
update_muta <- function() {

  # Update RR for each SNP
  for (idx in 1:snp.nums) {

    # Initialize rr.mh.chain with length = (rr.mh.updates+1)
    mh.chain <- c(mr.current[idx], rep(NA,mr.mh.updates))
    # store acceptance results
    mr.mh.accept <- rep(NA, mr.mh.updates)

    # parameter
    hh <-  HH[iter+1, idx]

    for (mh.idx in 1:mr.mh.updates ) {

      # current MH state
      muti <- mh.chain[mh.idx]

      # propose new value
      shape.curr <- ( mr_s1 *muti+ 1) / (1-muti)
      # rate.curr_mr <- mr_s1
      new.muti <- rbeta(1, shape.curr, mr_s1)

      # rate 1
      tmp1 <- multinom_ratio( gam = ri.current[idx],
                              alpha = af.current[idx],
                              miu = new.muti,
                              gmodel = gmdl)   /
        multinom_ratio( gam = ri.current[idx],
                        alpha = af.current[idx],
                        miu = muti,
                        gmodel = gmdl)
      sumlog <- sum(ssample[idx,]* log(tmp1))
      rate1 <- exp(sumlog)


      # rate 2 - prior distribution
      rate2 <- ifelse(hh,
                      dbeta(new.muti, shape1 = MR_alpha_1[idx], shape2 = MR_beta_1[idx])/
                        dbeta(muti, shape1 = MR_alpha_1[idx], shape2 = MR_beta_1[idx])
                      ,
                      dbeta(new.muti, shape1 = MR_alpha_0[idx], shape2 = MR_beta_0[idx])/
                        dbeta(muti, shape1 = MR_alpha_0[idx], shape2 = MR_beta_0[idx]) )


      # rate 3 - proposal jump probability
      shape.new_mr <- ( mr_s1 *new.muti+ 1) / (1-new.muti)

      rate3 <- dbeta(muti, shape.new_mr, mr_s1) /
        dbeta(new.muti, shape.curr, mr_s1)

      rate <- rate1 * rate2 * rate3
      accept <- rbinom(1,1,min(1,rate))

      # update next MH state
      mh.chain[mh.idx+1] <- new.muti*accept + muti*(1-accept)
      mr.mh.accept[mh.idx] <- accept

    }

    # update next Gibbs
    mr.current[idx] <<- mean(tail(mh.chain, 5))
    mmutate[iter+1, idx] <<- mr.current[idx]
    # average acceptance rate
    mr_accept_rt[iter, idx] <<- mean(mr.mh.accept)

  }

  return()
}
