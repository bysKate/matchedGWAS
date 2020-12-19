


# Data --------------------------------------------------------------------

assign(
  x = 'data.maf',
  value = ((ssample[,4]+ssample[,5]+ssample[,6])/2 + ssample[,7]+ssample[,8]+ssample[,9])/sum(ssample[1,]),
  envir = .GlobalEnv
)

assign(
  x = 'data.mr',
  value = (ssample[,2] + ssample[,3] + ssample[,4] + ssample[,6] + ssample[,7] + ssample[,8]) / rowSums(ssample) /2,
  envir = .GlobalEnv
)

# If there are no MR observed, adjust for MR_obs=0
data.mr[which(data.mr==0)] <- 0.003

# Prob b & Gene seg status G ----------------------------------------------------------

assign(x = 'alpha_b', value = 1.1, envir = .GlobalEnv)
assign(x = 'beta_b', value = 2.2, envir = .GlobalEnv)


# SNP status -----------------------------------------------------------------

if (dim(ssample)[1] > 1) {

  # Multiple marker analysis

  # P(H=1 | G=1)
  assign(x = 'p1', value = 1/log10(snp.size + 9), envir = .GlobalEnv)
  # P(H=1 | G=0)
  assign(x = 'p0', value = 0.3, envir = .GlobalEnv)

} else {

  # Single marker analysis

  # P(H=1 | G=1)
  assign(x = 'p1', value = 0.99, envir = .GlobalEnv)
  # P(H=1 | G=0)
  assign(x = 'p0', value = 0.01, envir = .GlobalEnv)

}


# Relative Risk --------------------------------------------------------------------

# Metropolis-Hasting iteration
assign(x = 'rr.mh.updates', value = 39, envir = .GlobalEnv)

# R | H=0 ~ Gamma(RR_shape_0, RR_rate_0)
assign(x = 'RR_shape_0', value = rep(3, snp.size), envir = .GlobalEnv)
assign(x = 'RR_rate_0', value = rep(3, snp.size), envir = .GlobalEnv)

# R | H=1 ~ Gamma(RR_shape_1, RR_rate_1)
assign(x = 'RR_shape_1', value = rep(8, snp.size), envir = .GlobalEnv)
assign(x = 'RR_rate_1', value = rep(3, snp.size), envir = .GlobalEnv)

#  proposal coefficient -- Gamma(shape=1+shape.coef * risk, rate=rate.coef)
assign(x = 'shape.coef', value = 3, envir = .GlobalEnv)
assign(x = 'rate.coef', value = 3, envir = .GlobalEnv)




# Allele Frequency ------------------------------------------------------------------

# Metropolis-Hasting iteration
assign(x = 'af.mh.updates', value = 29, envir = .GlobalEnv)

# MAF | H=0 ~ Beta(MAF_alpha_0, MAF_beta_0)
n1 <- rep(10, snp.size)
assign(x = 'MAF_alpha_0', value = n1 + 1/3, envir = .GlobalEnv)
assign(x = 'MAF_beta_0', value = n1/data.maf - MAF_alpha_0 + 2/3, envir = .GlobalEnv)


# MAF | H=1 ~ Beta(MAF_alpha_1, MAF_beta_1)
mr0 <- 0.7 * data.mr
af1.center <- (3*data.maf*(2*mr0 - 1) - ((3*data.maf*(2*mr0 - 1) + 7/2)^2 - 12*data.maf*(2*mr0 - 1)*((3*mr0)/2 + 1))^(1/2) + 7/2)/(3*(2*mr0 - 1))

n2 <- rep(4, snp.size)

assign(x = 'MAF_alpha_1', value = n2 + 1/3, envir = .GlobalEnv)
assign(x = 'MAF_beta_1', value = n2/af1.center - MAF_alpha_1 + 2/3, envir = .GlobalEnv)


# proposal coefficient -- shape.curr=( af_coef *afi + 1) / (1-afi)
assign(x = 'MAF_shape_1', value = 120, envir = .GlobalEnv)




# Mutation Rate -----------------------------------------------------------

# Metropolis-Hasting iteration
assign(x = 'mr.mh.updates', value = 19, envir = .GlobalEnv)

#  MR | H = 0 ~ Beta(MR_alpha_0, MR_beta_0)
assign(
  x = 'MR_alpha_0',
  value = rep(25 + 1/3,  snp.size),
  envir = .GlobalEnv
)

assign(
  x = 'MR_beta_0',
  value = MR_alpha_0 / data.mr - MR_alpha_0 + 2/3,
  envir = .GlobalEnv
)


# MR | H = 1 ~ Beta(MR_alpha_1, MR_beta_1)
mr1.center <- (3* data.mr *(2* af1.center  - 1) - ((3* data.mr *(2* af1.center  - 1) + 7/2)^2 - 12* data.mr *(2* af1.center  - 1)*((3* af1.center )/2 + 1))^(1/2) + 7/2)/(3*(2* af1.center  - 1))

assign(
  x = 'MR_alpha_1',
  value = rep(10 + 1/3,  snp.size),
  envir = .GlobalEnv
)

assign(
  x = 'MR_beta_1',
  value =MR_alpha_1 / mr1.center - MR_alpha_1 + 2/3,
  envir = .GlobalEnv
)


# proposal coefficient
assign(x = 'mr_s1', value = 1000, envir = .GlobalEnv)





