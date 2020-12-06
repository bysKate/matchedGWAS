#' Update SNP status parameter H.
#' @description For J SNPs on the gene segment, iteratively update SNP status (H_1, H_2, ...., H_J).
#' Each SNP status is updated based on the most current state of all other parameters.
#'
#' @return
#' @export
#'
#' @examples
update_H <- function() {
  # parameters
  part1 <- p1^gg.current * p0^(1-gg.current)
  part0 <- (1-p1)^gg.current * (1-p0)^(1-gg.current)

  # one-by-one
  for (indx in 1:snp.nums) {
    # parameters
    gam1 <- dgamma(ri.current[indx], RR_shape_1[indx], RR_rate_1[indx])   # R | H=1
    gam0 <- dgamma(ri.current[indx], RR_shape_0[indx], RR_rate_0[indx])   # R | H=0

    maf1 <- dbeta(af.current[indx], MAF_alpha_1[indx], MAF_beta_1[indx]) # AF| H=1
    maf0 <- dbeta(af.current[indx], MAF_alpha_0[indx], MAF_beta_0[indx]) # AF| H=0

    mut1 <- dbeta(mr.current[indx], MR_alpha_1[indx], MR_beta_1[indx] )  # MR| H=1
    mut0 <- dbeta(mr.current[indx], MR_alpha_0[indx], MR_beta_0[indx] )  # MR| H=0



    # prob
    prob_hh <- ( part1 * gam1 * maf1 * mut1) / (
      part0 * gam0 * maf0 * mut0 + part1 * gam1 * maf1 * mut1  )

    # next state
    new.H <- rbinom(1,1,prob = prob_hh)
    HH[iter+1,indx] <<- new.H
  }

  assign('hh.current', HH[iter+1,], envir = .GlobalEnv)

  return()
}
