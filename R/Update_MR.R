#' Update mutation rate parameter.
#'
#' @return
#' @export
#'
#' @examples
update_muta <- function() {

  # by SNP
  for (idx in 1:snp.nums) {

    mh.chain <- c(mr.current[idx], rep(NA,mr.mh.updates))
    # parameter
    hh <-  HH[iter+1, idx]

    for (mh.idx in 1:mr.mh.updates ) {
      muti <- mh.chain[mh.idx]

      # propose new value
      # shape.curr <- ( mr_s1 *muti+ 1) / (1-muti)
      rate.curr_mr <- mr_s1 * (1-muti) / muti - 1/3/muti + 2/3
      new.muti <- rbeta(1, mr_s1, rate.curr_mr)

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
      rate.new_mr <- mr_s1 * (1-new.muti) / new.muti - 1/3/new.muti + 2/3

      rate3 <- dbeta(muti, mr_s1, rate.new_mr) /
        dbeta(new.muti, mr_s1, rate.curr_mr)

      rate <- rate1 * rate2 * rate3
      accept <- rbinom(1,1,min(1,rate))

      mh.chain[mh.idx+1] <- new.muti*accept + muti*(1-accept)
    }
    # ts.plot(mh.chain)
    # update
    mr.current[idx] <<- mh.chain[mr.mh.updates + 1]
    mmutate[iter+1, idx] <<- mr.current[idx]
  }

  return()
}
