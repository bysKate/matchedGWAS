#' Update allele frequency parameter.
#'
#' @return
#' @export
#'
#' @examples
update_al <- function() {

  # by SNP
  for (idx in 1:snp.nums) {

    # Initialize af.mh.chain with length = af.mh.updates + 1
    af.mh.chain <- c(af.current[idx], rep(NA,af.mh.updates))

    # parameters
    hh <-  HH[iter+1, idx]

    for (mh.idx in 1:af.mh.updates ) {
      ali <-  af.mh.chain[mh.idx]

      # propose rate
      rate.curr <- MAF_shape_1 * (1-ali) / ali - 1/3/ali + 2/3

      new.ali <- rbeta(1, MAF_shape_1, rate.curr)

      # rate 1
      tmp1 <- multinom_ratio(gam = ri.current[idx],
                             alpha = new.ali,
                             miu = mr.current[idx],
                             gmodel = gmdl)   /
        multinom_ratio(gam = ri.current[idx],
                       alpha = ali,
                       miu = mr.current[idx],
                       gmodel = gmdl)

      sumlog <- sum(ssample[idx,]* log(tmp1))
      rate1 <- exp(sumlog)

      # rate 2 - prior
      rate2 <- ifelse(hh,
                      dbeta(new.ali, shape1 = MAF_alpha_1[idx], shape2 = MAF_beta_1[idx])/
                        dbeta(ali, shape1 = MAF_alpha_1[idx], shape2 = MAF_beta_1[idx])
                      ,
                      dbeta(new.ali, shape1 = MAF_alpha_0[idx], shape2 = MAF_beta_0[idx])/
                        dbeta(ali, shape1 = MAF_alpha_0[idx], shape2 = MAF_beta_0[idx]) )

      # rate 3 - proposal jump probability
      rate.new <- MAF_shape_1 * (1-new.ali) / new.ali - 1/3/new.ali + 2/3

      rate3 <- dbeta( ali, MAF_shape_1, rate.new) /
        dbeta( new.ali, MAF_shape_1, rate.curr)

      # accept ratio
      rate <- rate1 * rate2 * rate3
      accept <- rbinom(1,1,min(1,rate))

      af.mh.chain[mh.idx+1] <- new.ali * accept + ali * (1-accept)

    }
    # ts.plot(af.mh.chain)
    # update
    af.current[idx] <<- af.mh.chain[af.mh.updates+1]
    aallele[iter+1, idx] <<- af.current[idx]
  }

  return()
}
