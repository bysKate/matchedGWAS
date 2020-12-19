#' Update allele frequency parameter.
#'
#' @description For J SNPs on the segment, iteratively update relative risk parameter (A_1, ..., A_J). The next state is updated based on the most current state of other parameters.
#' @return
#' @export
#'
#' @examples
update_al <- function() {

  # Update Allele freq for each SNP
  for (idx in 1:snp.nums) {

    # Initialize af.mh.chain with length = af.mh.updates + 1
    af.mh.chain <- c(af.current[idx], rep(NA,af.mh.updates))
    # store acceptance results
    af.mh.accept <- rep(NA, af.mh.updates)

    # parameters
    hh <-  HH[iter+1, idx]

    for (mh.idx in 1:af.mh.updates ) {
      ali <-  af.mh.chain[mh.idx]

      # propose rate
      rate.curr <-  MAF_shape_1 * (1-ali) / ali

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
      rate.new <-  MAF_shape_1 * (1-new.ali) / new.ali

      rate3 <- dbeta( ali, MAF_shape_1, rate.new) /
        dbeta( new.ali, MAF_shape_1, rate.curr)

      # accept ratio
      rate <- rate1 * rate2 * rate3
      accept <- rbinom(1,1,min(1,rate))

      # update next MH state
      af.mh.chain[mh.idx+1] <- new.ali * accept + ali * (1-accept)
      af.mh.accept[mh.idx] <- accept

    }

    # update
    af.current[idx] <<- mean(tail(af.mh.chain, 10))
    aallele[iter+1, idx] <<- af.current[idx]
    # average acceptance rate
    af_accept_rt[iter, idx] <<- mean(af.mh.accept)

  }

  return()
}
