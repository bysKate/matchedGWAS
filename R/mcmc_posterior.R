#' Markov Chain Monte Carlo (MCMC) simulations.
#'
#' @param cur_gid
#' @param saveResult
#' @param swap
#' @param swap_cutoff
#'
#' @return
#' @export
#'
#' @examples
mcmc_posterior <- function(
  cur_gid,
  swap = T,
  swap_cutoff = 0.8,
  saveResult = T
) {

  #  load  dataset
  ssample <<- as.matrix(
    allsnp_list %>% filter(gseg_id==cur_gid) %>% select(X00:X22),
    ncol = 9
  )

  # reverse risk allele
  if (isTRUE(swap)) {
    rec_idx <- (allsnp_list$rr_mle[allsnp_list$gseg_id==cur_gid] < swap_cutoff)
    ssample[rec_idx,] <<- t(
      apply(
        matrix(ssample[rec_idx,], byrow = T, ncol = 9), 1, rev
      )
    )
  }

  gmdl <<- 'add0'
  snp.size <<- dim(ssample)[1]


  # source codes
  setwd(curr_dir)
  source('params.R')


  # mcmc 1
  update_mcmc(MCchains_len, use_seed = 1)
  mclist1 <- esti.list; rm(esti.list, envir = .GlobalEnv)

  # mcmc 2
  update_mcmc(MCchains_len, use_seed = 2)
  mclist2 <- esti.list; rm(esti.list, envir = .GlobalEnv)

  # mcmc 3
  update_mcmc(MCchains_len, use_seed = 3)
  mclist3 <- esti.list; rm(esti.list, envir = .GlobalEnv)

  # combine results to list
  mcmcresult <- list('mclist1' = mclist1,
                     'mclist2' = mclist2,
                     'mclist3' = mclist3
  )



  if (isTRUE(saveResult)) {

    # tmp directory
    setwd(result_path)

    if (!dir.exists(result_folder)) { dir.create(result_folder) }
    setwd(result_folder)

    if (!dir.exists('results')) { dir.create('results') }
    setwd('results')

    # save
    saveRDS(mcmcresult, file = paste('new_Gid_',cur_gid,'.rds',sep = ''))

    return()
  } else {
    return(mcmcresult)
  }

}
