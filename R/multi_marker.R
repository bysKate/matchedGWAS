#' Apply multiple-marker analysis using hierarchical Bayesian model.
#'
#' @param cur_gid
#' @param saveResult
#'
#' @return
#' @export
#'
#' @examples
multi_MCMC_data <-function(cur_gid, saveResult = T) {
  # Notes


  #  load  dataset
  ssample <<- as.matrix(allsnp_list[[cur_gid]][,4:12], ncol = 9)

  # additional Rev
  rev_ind <- allsnp_list[[cur_gid]]$snpid %in% rev_SNP

  # # # #

  gmdl <<- 'add0'
  snp.size <<- dim(ssample)[1]

  if (sum(rev_ind) > 1) {
    ssample[rev_ind,] <<- t(apply(ssample[rev_ind,], 1, rev))
  } else if (sum(rev_ind) == 1) {
    ssample[rev_ind,] <<- rev(ssample[rev_ind,])
  }

  # source codes
  setwd(curr_dir)
  source('params_v4.R')


  # mcmc 1
  update_mcmc(MCchains_len,initl = 1)
  mclist1 <- esti.list; rm(esti.list, envir = .GlobalEnv)

  # mcmc 2
  update_mcmc(MCchains_len,initl = 2)
  mclist2 <- esti.list; rm(esti.list, envir = .GlobalEnv)

  # mcmc 3
  update_mcmc(MCchains_len,initl = 3)
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
