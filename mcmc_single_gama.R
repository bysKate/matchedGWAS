##############   mcmc update (Single SNP) ####################
update_single <- function(sim.steps,initl=NA) {
  # esti.list is a list of estimation objects
  # ssample: a snp.nums * 9 matrix
  # snp.nums: number of SNPs for analysis
  
  snp.nums <<- dim(ssample)[1]
  rrisk <<- matrix(ncol = snp.nums, nrow = sim.steps)
  aallele <<- matrix(ncol = snp.nums, nrow = sim.steps)
  mmutate <<- matrix(ncol = snp.nums, nrow = sim.steps)
  
  for (snp.id in 1:snp.nums) {
    idx <<- snp.id
    countwhl <<- 0
    
    
    while ( any(is.na(rrisk[20,idx]), (aallele[20,idx]==aallele[19,idx])) && (countwhl<50)) {
      countwhl <<- countwhl+1
      #print(paste('count',countwhl))
      #initialize
      if (is.na(initl)) {
        # initialisation
        aallele[1,idx] <<- runif(1,0,1)  #need global assignments
        mmutate[1,idx] <<- runif(1,min = 0,max = 0.1)
        rrisk[1,idx] <<- rgamma(1,shape = 5,rate = 5)
      } else {
        aallele[1,idx] <<- 0.35
        mmutate[1,idx] <<- 0.0007
        rrisk[1,idx] <<- 1
      }
      
      for ( iter.local in 1:20 ) {
        iter <<- iter.local
        ri.current <<- rrisk[iter, idx]
        af.current <<- aallele[iter, idx]
        mr.current <<- mmutate[iter, idx]
        # Metropolis-Hasting, randomize update orders
        ordr=sample(x = 3,size = 3)
        for (oo in ordr) {
          if (oo==1) { update_risk_single() } else if (oo==2) { update_al_single() } else { update_muta_single() }
        }
      }
    }
    
    
    if (countwhl < 50) {
      for ( iter.local in 21:(sim.steps-1) ) {
        iter <<- iter.local
        ri.current <<- rrisk[iter, idx]
        af.current <<- aallele[iter, idx]
        mr.current <<- mmutate[iter, idx]
        # Metropolis-Hasting, randomize update orders
        ordr=sample(x = 3,size = 3)
        for (oo in ordr) {
          if (oo==1) { update_risk_single() } else if (oo==2) { update_al_single() } else { update_muta_single() }
        }
      }
    } else { print('count >= 50, update stopped!!!') }
  }
  
  
  # combine list
  esti.list <<- list('rrisk'=rrisk, 'aallele'=aallele, 'mmutate'=mmutate)
}


update_risk_single <- function() {
  mh.chain=c(ri.current, rep(0,99))

  for (mh.idx in 1:99 ) {
    ri = mh.chain[mh.idx]
    shapenew=1 + 5 * ri
    new.ri=rgamma(1,shape = shapenew, rate = 5)
    #print(paste('risk',ri,'and new risk',new.ri))     
    
    tmp1 = multinom_ratio(miu = mr.current, alpha = af.current, gam = new.ri, gmodel = gmdl)/
      multinom_ratio(miu = mr.current, alpha = af.current, gam = ri, gmodel = gmdl)
    rate1 = prod(tmp1 ^ ssample[idx,])   # product of 9 cells
    #print((paste('rate1',rate1)))
    if(is.na(rate1)) {mh.chain[mh.idx+1]=ri; next;}
    
    # rate 3
    rate3 = dgamma(ri,shape = 1 + 5 * new.ri, rate = 5) / dgamma(new.ri,shape = shapenew, rate = 5)
    #print(paste('rate3',rate3))
    rate = rate1 * rate3
    accept = rbinom(1,1,min(1,rate))
    #print(paste("-------------- Accept rate for 'Risk",idx,"' is ",rate, " !", sep = ''))
    
    mh.chain[mh.idx+1] = new.ri * accept + ri * (1-accept)
  }
  
  # ts.plot(mh.chain)
  ri.current <<- mh.chain[sample(90:100,1)]
  rrisk[iter+1, idx] <<- ri.current
}



update_al_single <- function() {
  mh.chain=c(af.current, rep(0,99))
  
  for (mh.idx in 1:99 ) {
    ali = mh.chain[mh.idx]
    shapenew=( 100*ali + 1) / (1-ali)
    new.ali=rbeta(1,shapenew,102)   # let shape2 be 102 
    # rate 1
    tmp1 = multinom_ratio(miu = mr.current, gam = ri.current, alpha = new.ali, gmodel = gmdl)/
      multinom_ratio(miu = mr.current, gam = ri.current, alpha = ali, gmodel = gmdl)
    rate1=prod(tmp1 ^ ssample[idx,])
    
    if(is.na(rate1)) {mh.chain[mh.idx+1]=ali; next}
    
    # rate 3
    rate3=dbeta(new.ali, shapenew, 102) / dbeta(ali, ( 100*new.ali + 1)/(1-new.ali), 102)
    
    rate=rate1 / rate3
    accept=rbinom(1,1,min(1,rate))
    #print(paste("-------------- Accept rate for 'Allele",idx,"' is ",rate, " !", sep = ''))
    
    mh.chain[mh.idx+1] = new.ali*accept + ali*(1-accept)
  }
  
  # ts.plot(mh.chain)
  af.current <<- mh.chain[sample(90:100,1)]
  aallele[iter+1, idx] <<- af.current
}



update_muta_single <- function() {
  mh.chain=c(mr.current, rep(0,99))
  
  for (mh.idx in 1:99 ) {
    muti = mh.chain[mh.idx]
    shapenew=( 1000*muti+ 1) / (1-muti)    
    new.muti=rbeta(1,shapenew,1002) 		# let shape2=1002
    # rate 1
    tmp1 = multinom_ratio( alpha = af.current, gam = ri.current, miu = new.muti, gmodel = gmdl)/
      multinom_ratio( alpha = af.current, gam = ri.current, miu = muti, gmodel = gmdl)
    rate1 = prod(tmp1 ^ ssample[idx,])
    if(is.na(rate1)) {mh.chain[mh.idx+1]=muti; next}
    #print(paste('rate1',rate1))
    
    # rate 3
    rate3 = dbeta(new.muti, shapenew, 1002) / dbeta(muti, ( 1000*new.muti+ 1)/(1-new.muti), 1002)
    #print(paste('rate3',rate3))
    rate=rate1 / rate3
    accept=rbinom(1,1,min(1,rate))
    #print(paste("-------------- Accept rate for 'Mutation",idx,"' is ",rate, " !", sep = ''))
    
    mh.chain[mh.idx+1] = new.muti*accept + muti*(1-accept)
  }
  # ts.plot(mh.chain)
  mr.current <<- mh.chain[sample(90:100,1)]
  mmutate[iter+1, idx] <<- mr.current
}