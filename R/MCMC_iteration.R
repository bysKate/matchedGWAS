#' MCMC iteration.
#'
#' @param itt iteration numbers
#' @param initl initialization number
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
update_mcmc <- function(itt, initl=NA) {

  assign('sim.steps',
         itt,
         envir = .GlobalEnv)
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


  # initialisation
  if (initl==1) {
    set.seed(1)

    bb[1] <<- 0.3
    ggstatus[1] <<- 0
    HH[1,] <<- 0

    rrisk[1,] <<- 1
    aallele[1,] <<- data.maf
    mmutate[1,] <<- 0.005

  } else if (initl==2){
    set.seed(2)

    bb[1] <<- 0.5
    ggstatus[1] <<- 1
    HH[1,] <<- 1

    rrisk[1,] <<- 3         # large risk
    aallele[1,] <<- af1.center
    mmutate[1,] <<- 0.005

  } else {
    set.seed(3)

    bb[1] <<- 0.3
    ggstatus[1] <<- 0
    HH[1,] <<- 1

    rrisk[1,] <<- 1.8       # moderate risk
    aallele[1,] <<- (af1.center+data.maf)/2
    mmutate[1,] <<- 0.005

  }

  # iterations
  for ( iter.local in 1:(sim.steps-1) ) {

    # assign current value
    assign('iter', iter.local, envir = .GlobalEnv)
    assign('bb.current', bb[iter.local] , envir = .GlobalEnv)
    assign('gg.current', ggstatus[iter.local], envir = .GlobalEnv)
    assign('hh.current', HH[iter.local,], envir = .GlobalEnv)
    assign('ri.current', rrisk[iter.local, ], envir = .GlobalEnv)
    assign('af.current', aallele[iter.local, ], envir = .GlobalEnv)
    assign('mr.current', mmutate[iter.local, ], envir = .GlobalEnv)

    if (is.na(bb.current) & is.na(gg.current)) {
      esti.list <<- list('bb'=bb,
                         'ggstatus'=ggstatus,
                         'HH'=HH,
                         'rrisk'=rrisk,
                         'aallele'=aallele,
                         'mmutate'=mmutate
      )
      return()
    }
    # Gibbs update
    update_b1()
    update_gstatus()
    update_H()

    # Metropolis-Hasting, randomize update orders
    ordr=sample(x = 3,size = 3)
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
  esti.list <<- list('bb'=bb,
                     'ggstatus'=ggstatus,
                     'HH'=HH,
                     'rrisk'=rrisk,
                     'aallele'=aallele,
                     'mmutate'=mmutate
  )

}
