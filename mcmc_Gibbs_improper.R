multinom_ratio <- function(miu, alpha, gam, gmodel) {
  # gamma distribution(shape, rate)
  # ssample: snp.nums * 9 matrix 
  
  if (gmodel=='add0') {
    
    increm=(gam-1)/2
    q00 = (1-alpha)*(1-miu)*(1-alpha)*(1-miu)
    q01 = (1-alpha)*(1-alpha)*2*miu*(1-miu)* (1+increm) # if R1=gam, then R2 could be negative!
    q02 = (1-alpha)*(1-alpha)*miu*miu* gam
    
    q10 = 2*(1-alpha)*(1-miu)*alpha*miu
    q11 = 2*alpha*(1-alpha)*(2*miu*miu-2*miu+1)* (1+increm)
    q12 = 2*alpha*(1-alpha)*miu*(1-miu)* gam
    
    q20 = alpha*alpha*miu*miu
    q21 = alpha*alpha*2*miu*(1-miu)* (1+increm)
    q22 = alpha*alpha*(1-miu)*(1-miu)* gam
    
  } else if (gmodel=='rec') {
    
    q00 = (1-alpha)*(1-miu)*(1-alpha)*(1-miu)
    q01 = (1-alpha)*(1-alpha)*2*miu*(1-miu)*1
    q02 = (1-alpha)*(1-alpha)*miu*miu*gam 
    
    q10 = 2*(1-alpha)*(1-miu)*alpha*miu
    q11 = 2*alpha*(1-alpha)*(2*miu*miu-2*miu+1)*1
    q12 = 2*alpha*(1-alpha)*miu*(1-miu)*gam 
    
    q20 = alpha*alpha*miu*miu
    q21 = alpha*alpha*2*miu*(1-miu)* 1
    q22 = alpha*alpha*(1-miu)*(1-miu)*gam 
    
  } else if (gmodel=='dom') {
    
    q00 = (1-alpha)*(1-miu)*(1-alpha)*(1-miu)
    q01 = (1-alpha)*(1-alpha)*2*miu*(1-miu)*gam
    q02 = (1-alpha)*(1-alpha)*miu*miu*gam 
    
    q10 = 2*(1-alpha)*(1-miu)*alpha*miu
    q11 = 2*alpha*(1-alpha)*(2*miu*miu-2*miu+1)*gam
    q12 = 2*alpha*(1-alpha)*miu*(1-miu)*gam 
    
    q20 = alpha*alpha*miu*miu
    q21 = alpha*alpha*2*miu*(1-miu)*gam
    q22 = alpha*alpha*(1-miu)*(1-miu)*gam 
    
  } else if (gmodel=='mul') {
    
    q00 = (1-alpha)*(1-miu)*(1-alpha)*(1-miu)
    q01 = (1-alpha)*(1-alpha)*2*miu*(1-miu)* gam
    q02 = (1-alpha)*(1-alpha)*miu*miu* gam^2
    
    q10 = 2*(1-alpha)*(1-miu)*alpha*miu
    q11 = 2*alpha*(1-alpha)*(2*miu*miu-2*miu+1)* gam
    q12 = 2*alpha*(1-alpha)*miu*(1-miu)* gam^2
    
    q20 = alpha*alpha*miu*miu
    q21 = alpha*alpha*2*miu*(1-miu)* gam
    q22 = alpha*alpha*(1-miu)*(1-miu)* gam^2
    
  } else if (gmodel=='dou') {
    q00 = (1-alpha)*(1-miu)*(1-alpha)*(1-miu)
    q01 = (1-alpha)*(1-alpha)*2*miu*(1-miu)*gam
    q02 = (1-alpha)*(1-alpha)*miu*miu*gam *2 
    
    q10 = 2*(1-alpha)*(1-miu)*alpha*miu
    q11 = 2*alpha*(1-alpha)*(2*miu*miu-2*miu+1)*gam
    q12 = 2*alpha*(1-alpha)*miu*(1-miu)*gam *2
    
    q20 = alpha*alpha*miu*miu
    q21 = alpha*alpha*2*miu*(1-miu)*gam
    q22 = alpha*alpha*(1-miu)*(1-miu)*gam *2
  }
  
  deno =q00+q01+q02+q10+q11+q12+q20+q21+q22
  proportion = rbind(q00/deno,q01/deno,q02/deno,
                     q10/deno,q11/deno,q12/deno,
                     q20/deno,q21/deno,q22/deno)
  #print(deno)
  return(proportion)
  
} 




########################## Metropolis Hasting algorithm #########################
# At step of 'iter'
update_b1 <- function() {
  shap1 = ggstatus[iter] +1   #alpha_b = 1
  shap2 = 2-ggstatus[iter]   #beta_b = 1
  bb[iter+1] <<- rbeta(1,shape1 = shap1,shape2 = shap2)
}

update_gstatus <- function() { 
  
  # parameters
  b=bb[iter+1]
  gstatus=ggstatus[iter]
  sum.H=sum(HH[iter,])     # sum of 'iter' row of HH
  numer1 = b * ((1/2)/(1-1/2))^sum.H * (1-1/2)^snp.nums
  # gibbs
  prob_G = numer1/ ( numer1 + ( (1-b) * (1/10/(1-1/10))^sum.H * (1-1/10)^snp.nums ) )
  ggstatus[iter+1] <<- rbinom(1,1,prob = prob_G)
}

update_H <- function() {
  # parameters
  gstatus = ggstatus[iter+1]
  part1 = (1/2)^gstatus * 1/10^(1-gstatus)
  part0 = (1-1/2)^gstatus * (1-1/10)^(1-gstatus)
  
  # one-by-one
  for (indx in 1:snp.nums) {
    # parameters
    gam1 = 1   # improper
    gam0 = 1   # improper
    
    prob_hh = (gam1*part1) / (gam0*part0+gam1*part1)
    # gibbs
    HH[iter+1,indx] <<- rbinom(1,1,prob = prob_hh)
  }
}



update_risk <- function() {
  
  # one-by-one
  for (idx in 1:snp.nums) {
    
    hh = HH[iter+1, idx]
    
    mh.chain=c(ri.current[idx], rep(0,99))
    for (mh.idx in 1:99 ) {
      ri = mh.chain[mh.idx]
      shapenew=1 + 5 * ri
      new.ri=rgamma(1,shape = shapenew, rate = 5)
      #print(paste('risk',ri,'and new risk',new.ri))     
      tmp1 = multinom_ratio(miu = mr.current[idx], 
                            alpha = af.current[idx], 
                            gam = new.ri, 
                            gmodel = gmdl)   /
        multinom_ratio(miu = mr.current[idx], 
                       alpha = af.current[idx], 
                       gam = ri, 
                       gmodel = gmdl)
      rate1 = prod(tmp1 ^ ssample[idx,])   # product of 9 cells
      #print((paste('rate1',rate1)))
      if(is.na(rate1)) {mh.chain[mh.idx+1]=ri; next;}
      
      rate2 = 1   # improper
      
      rate3 = dgamma(ri,shape = 1 + 5 * new.ri, rate = 5) / dgamma(new.ri,shape = shapenew, rate = 5)
      #print(paste('rate3',rate3))
      rate = rate1 * rate2 * rate3
      accept = rbinom(1,1,min(1,rate))
      #print(paste("-------------- Accept rate for 'Risk",idx,"' is ",rate, " !", sep = ''))
      
      mh.chain[mh.idx+1] = new.ri * accept + ri * (1-accept)
    }
    
    # ts.plot(mh.chain)
    # sample & update
    ri.current[idx] <<- mh.chain[sample(90:100,1)]
    rrisk[iter+1, idx] <<- ri.current[idx]
  }
}





update_al <- function() {
  
  # one-by-one SNPs
  for (idx in 1:snp.nums) {
    mh.chain=c(af.current[idx], rep(0,99))
    
    for (mh.idx in 1:99 ) {
      ali = mh.chain[mh.idx]
      shapenew=( 100*ali/crng + 1) / (1-ali/crng)
      new.xi = rbeta(1,shapenew,102)
      new.ali = new.xi * crng   # let shape2 be 102 
      # rate 1
      tmp1 = multinom_ratio(miu = mr.current[idx], 
                            gam = ri.current[idx], 
                            alpha = new.ali, 
                            gmodel = gmdl)   /
        multinom_ratio(miu = mr.current[idx], 
                       gam = ri.current[idx], 
                       alpha = ali, 
                       gmodel = gmdl)
      rate1=prod(tmp1 ^ ssample[idx,])
      
      if(is.na(rate1)) {mh.chain[mh.idx+1]=ali; next}
      
      #rate2 = dbeta(new.ali,4,8) / dbeta(ali,4,8)
      rate2=1
      # rate 3
      rate3=dbeta(new.xi, shapenew, 102) / dbeta(ali/crng, ( 100 * new.xi + 1)/(1 - new.xi), 102)
      
      rate=rate1 *rate2/ rate3
      accept=rbinom(1,1,min(1,rate))
      #print(paste("-------------- Accept rate for 'Allele",idx,"' is ",rate, " !", sep = ''))
      
      mh.chain[mh.idx+1] = new.ali*accept + ali*(1-accept)
    }
    # ts.plot(mh.chain)
    # sample & update
    af.current[idx] <<- mh.chain[sample(90:100,1)]
    aallele[iter+1, idx] <<- af.current[idx]
  }
}



update_muta <- function() {
  
  # one-by-one SNPs
  for (idx in 1:snp.nums) {
    
    mh.chain=c(mr.current[idx], rep(0,99))
    
    for (mh.idx in 1:99 ) {
      muti = mh.chain[mh.idx]
      shapenew=( 1000*muti+ 1) / (1-muti)    
      new.muti=rbeta(1,shapenew,1002) 		# let shape2=1002
      # rate 1
      tmp1 = multinom_ratio( alpha = af.current[idx], 
                             gam = ri.current[idx], 
                             miu = new.muti, 
                             gmodel = gmdl)   /
        multinom_ratio( alpha = af.current[idx], 
                        gam = ri.current[idx], 
                        miu = muti, 
                        gmodel = gmdl)
      rate1 = prod(tmp1 ^ ssample[idx,])
      if(is.na(rate1)) {mh.chain[mh.idx+1]=muti; next}
      #print(paste('rate1',rate1))
      rate2=1
      
      rate3 = dbeta(new.muti, shapenew, 1002) / dbeta(muti, ( 1000*new.muti+ 1)/(1-new.muti), 1002)
      #print(paste('rate3',rate3))
      rate=rate1 * rate2 / rate3
      accept=rbinom(1,1,min(1,rate))
      #print(paste("-------------- Accept rate for 'Mutation",idx,"' is ",rate, " !", sep = ''))
      
      mh.chain[mh.idx+1] = new.muti*accept + muti*(1-accept)
    }
    # ts.plot(mh.chain)
    # update
    mr.current[idx] <<- mh.chain[sample(90:100,1)]
    mmutate[iter+1, idx] <<- mr.current[idx]
  }
}




######################   Main Function  ##########################
update_mcmc <- function(sim.steps,initl=NA) {
  # esti.list is a list of estimation objects
  # ssample: a snp.nums * 9 matrix
  # snp.nums: number of SNPs for analysis
  
  snp.nums <<- dim(ssample)[1]
  aallele <<- matrix(ncol = snp.nums, nrow = sim.steps)
  mmutate <<- matrix(ncol = snp.nums, nrow = sim.steps )
  rrisk <<- matrix(ncol = snp.nums, nrow = sim.steps )
  HH <<- matrix(ncol = snp.nums, nrow = sim.steps )
  maf.di <<- ((ssample[,2]+ssample[,5]+ssample[,8])/2 + ssample[,3]+ssample[,6]+ssample[,9])/sum(ssample[1,])
  crng <<- 0.5
  # range will be between 0.5, if assume
  
  
  countwhl <<- 1
  while ( any(c(is.na(rrisk[20,]), (aallele[20,]==aallele[19,]))) && (countwhl<200)) {
    
    if (is.na(initl)) {
      # initialisation
      aallele[1,] <<- runif(snp.nums,min = 0,max = maf.di)
      mmutate[1,] <<- runif(snp.nums,min = 0,max = 0.1)
      rrisk[1,] <<- rgamma(snp.nums,shape = 5,rate = 5)
      HH[1,] <<- rbinom(snp.nums,1,0.5)
      ggstatus <<- c(rbinom(1,1,0.5) , rep(0, sim.steps-1))
      bb <<- c(runif(1) , rep(0, sim.steps-1))
    } else if (initl=='null'){
      aallele[1,] <<- runif(snp.nums,min = 0,max = maf.di)
      mmutate[1,] <<- 0.0007
      rrisk[1,] <<- 1
      HH[1,] <<- 0
      ggstatus <<- rep(0, sim.steps)
      bb <<- c(0.1 , rep(0, sim.steps-1))
    } else {
      aallele[1,] <<- runif(snp.nums,min = 0,max = maf.di)
      mmutate[1,] <<- 0.0002
      rrisk[1,] <<- 3.9
      HH[1,] <<- 1
      ggstatus <<- rep(1, sim.steps)
      bb <<- c(0.7 , rep(0, sim.steps-1))
      
    }
    
    for ( iter.local in 1:20 ) {
      iter <<- iter.local
      ri.current <<- rrisk[iter, ]
      af.current <<- aallele[iter, ]
      mr.current <<- mmutate[iter, ]
      
      # Gibbs update
      update_b1()
      update_gstatus()
      update_H()
      
      # Metropolis-Hasting, randomize update orders
      ordr=sample(x = 3,size = 3)
      for (oo in ordr) {
        if (oo==1) { update_risk() } else if (oo==2) { update_al() } else { update_muta() }
      }
      iter <<- iter+1
    }
    countwhl <<- countwhl+1
  }
  
  
  if (countwhl < 200) {
    for ( iter.local in 21:(sim.steps-1) ) {
      iter <<- iter.local
      ri.current <<- rrisk[iter, ]
      af.current <<- aallele[iter, ]
      mr.current <<- mmutate[iter, ]
      
      # Gibbs update
      update_b1()
      update_gstatus()
      update_H()
      
      # Metropolis-Hasting, randomize update orders
      ordr=sample(x = 3,size = 3)
      for (oo in ordr) {
        if (oo==1) { update_risk() } else if (oo==2) { update_al() } else { update_muta() }
      }
      iter <<- iter+1
    }
    
    esti.list <<- list('bb'=bb, 'ggstatus'=ggstatus, 'HH'=HH, 
                       'rrisk'=rrisk, 'aallele'=aallele, 'mmutate'=mmutate)
  } else {print('count > 2000, update stopped!!!')}
  
}




############## burning and thining ##################
# thining
# from acf plot, thining every 25 obs
thining <- function(sris,thinby=1) {
  if (is.vector(sris)) {
    sris.len=length(sris)
    tmp=seq(1,sris.len,by=thinby)
    sris.thin=sris[tmp]
  } else  if (is.matrix(sris)) {
    sris.len=dim(sris)[1]
    tmp=seq(1,sris.len,by=thinby)
    sris.thin=as.matrix(sris[tmp,])
  }
  return(list(sris.thin))
}

burning <- function(sris,percnt) {
  if (is.vector(sris)) {
    sris.len=length(sris)
    sris.start=round(sris.len*percnt)
    return(list(sris[sris.start:sris.len]))
  } else {
    sris.len=dim(sris)[1]
    sris.start=round(sris.len*percnt)
    return(list(as.matrix(sris[sris.start:sris.len,])))
  }
}




splitFacet <- function(x){
  facet_vars <- names(x$facet$params$facets)         # 1
  x$facet    <- ggplot2::ggplot()$facet              # 2
  datasets   <- split(x$data, x$data[facet_vars])    # 3
  new_plots  <- lapply(datasets,function(new_data) { # 4
    x$data <- new_data
    x})
}    

