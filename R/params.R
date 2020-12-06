


# Data --------------------------------------------------------------------

data.maf <<- ((ssample[,4]+ssample[,5]+ssample[,6])/2 + ssample[,7]+ssample[,8]+ssample[,9])/sum(ssample[1,])

data.mr <<- (ssample[,2] + ssample[,3] + ssample[,4] + ssample[,6] + ssample[,7] + ssample[,8]) / rowSums(ssample) /2

# If there are no MR observed, adjust for MR_obs=0
data.mr[which(data.mr==0)] <- 0.003

# Prob b & Gene seg status G ----------------------------------------------------------

alpha_b <<- 1.1
beta_b <<- 2.2


# SNP status -----------------------------------------------------------------

if (dim(ssample)[1] > 1) {

  # Multiple marker analysis
  p1 <<- 1/log10(snp.size + 9)    # P(H=1 | G=1)
  p0 <<- 0.3                      # P(H=1 | G=0)
} else {

  # Single marker analysis
  p1 <<- 1    # P(H=1 | G=1)
  p0 <<- 0    # P(H=1 | G=0)
}


# Relative Risk --------------------------------------------------------------------

# Metropolis-Hasting iteration
rr.mh.updates <-  79

# R | H=0 ~ Gamma(RR_shape_0, RR_rate_0)
RR_shape_0 <<- rep(3, snp.size)
RR_rate_0 <<- rep(3, snp.size)

# R | H=1 ~ Gamma(RR_shape_1, RR_rate_1)
RR_shape_1 <<- rep(8, snp.size)
RR_rate_1 <<- rep(3, snp.size)

#  proposal coefficient
shape.coef = 3
rate.coef = 3 # Gamma(shape=1+shape.coef * risk, rate=rate.coef)




# Allele Frequency ------------------------------------------------------------------

# Metropolis-Hasting iteration
af.mh.updates <-  49

# MAF | H=0 ~ Beta(MAF_alpha_0, MAF_beta_0)
n1 <- rep(10, snp.size)
MAF_alpha_0 <- n1 + 1/3
MAF_beta_0 <- n1/data.maf - MAF_alpha_0 + 2/3


# MAF | H=1 ~ Beta(MAF_alpha_1, MAF_beta_1)
mr0 <- 0.7 * data.mr
af1.center <- (3*data.maf*(2*mr0 - 1) - ((3*data.maf*(2*mr0 - 1) + 7/2)^2 - 12*data.maf*(2*mr0 - 1)*((3*mr0)/2 + 1))^(1/2) + 7/2)/(3*(2*mr0 - 1))

n2 <- rep(4, snp.size)

MAF_alpha_1 <- n2 + 1/3
MAF_beta_1 <- n2/af1.center - MAF_alpha_1 + 2/3


# proposal coefficient
MAF_shape_1 = 2  # shape.curr=( af_coef *afi + 1) / (1-afi)




# Mutation Rate -----------------------------------------------------------

# Metropolis-Hasting iteration
mr.mh.updates <- 29

#  MR | H = 0 ~ Beta(MR_alpha_0, MR_beta_0)
MR_alpha_0 <- rep(25 + 1/3,  snp.size)
MR_beta_0 <- MR_alpha_0 / data.mr - MR_alpha_0 + 2/3


# MR | H = 1 ~ Beta(MR_alpha_0, MR_beta_0)
mr1.center <- (3* data.mr *(2* af1.center  - 1) - ((3* data.mr *(2* af1.center  - 1) + 7/2)^2 - 12* data.mr *(2* af1.center  - 1)*((3* af1.center )/2 + 1))^(1/2) + 7/2)/(3*(2* af1.center  - 1))

MR_alpha_1 <- rep(10 + 1/3,  snp.size)
MR_beta_1 <- MR_alpha_1 / mr1.center - MR_alpha_1 +2/3


# proposal coefficient
mr_s1 <- 1000





