#' MCMC iteration.
#'
#' @param use_seed specify seed for initialization
#' @param itt iteration numbers
#'
#' @return esti.list:  a list of estimation objects
#' @export
#' @details
#' Global objects:
#' ssample: a matrix [snp.nums * 9]
#' snp.nums: number of SNPs for analysis
#' iter: current status (Gibbs sampler)
#'
#' @examples
update_mcmc <- function(itt, use_seed = 1) {


  assign('snp.nums',
         dim(ssample)[1],
         envir = .GlobalEnv)

  sapply(c('bb','ggstatus'),
         assign,
         rep(NA, itt),
         envir = .GlobalEnv)

  sapply(c('HH','rrisk','aallele','mmutate'),
         assign,
         matrix(ncol = snp.nums, nrow = itt),
         envir = .GlobalEnv)

  sapply(c('rr_accept_rt', 'af_accept_rt', 'mr_accept_rt'),
         assign,
         matrix(ncol = snp.nums, nrow = itt-1),
         envir = .GlobalEnv
  )


  # initialisation
  set.seed(use_seed)

  if (use_seed==1) {

    # no risk
    bb[1] <<- 0.2
    ggstatus[1] <<- 0
    HH[1,] <<- 0
    rrisk[1,] <<- 1
    aallele[1,] <<- data.maf
    mmutate[1,] <<- 0.005

  } else if (use_seed==2){

    # moderate risk
    bb[1] <<- 0.5
    ggstatus[1] <<- 0
    HH[1,] <<- 1
    rrisk[1,] <<- 2
    aallele[1,] <<- af1.center
    mmutate[1,] <<- 0.005

  } else {

    # large risk
    bb[1] <<- 0.7
    ggstatus[1] <<- 1
    HH[1,] <<- 1
    rrisk[1,] <<- 3
    aallele[1,] <<- af1.center * 0.5
    mmutate[1,] <<- 0.005

  }

  # Update current status
  assign('bb.current', bb[1] ,      envir = .GlobalEnv)
  assign('gg.current', ggstatus[1], envir = .GlobalEnv)
  assign('hh.current', HH[1,],      envir = .GlobalEnv)
  assign('ri.current', rrisk[1, ],  envir = .GlobalEnv)
  assign('af.current', aallele[1, ], envir = .GlobalEnv)
  assign('mr.current', mmutate[1, ], envir = .GlobalEnv)


  # iterations
  for ( iter.local in 1:(itt-1) ) {

    assign('iter', iter.local, envir = .GlobalEnv)

    # Gibbs update
    update_b1()
    update_gstatus()
    update_H()

    # Metropolis-Hasting, randomize update orders
    ordr <- sample(x = 3,size = 3)

    for (oo in ordr) {

      if (oo==1) {
        update_risk()
      } else if (oo==2) {
        update_al()
      } else {
        update_muta()
      }

    }

  }

  # combine list
  esti.list <<- list(
    'bb' = bb,
    'ggstatus' = ggstatus,
    'HH' = HH,
    'rrisk' = rrisk,
    'aallele' = aallele,
    'mmutate' = mmutate,
    'rr_accept_rt' = rr_accept_rt,
    'af_accept_rt' = af_accept_rt,
    'mr_accept_rt' = mr_accept_rt
  )

}
