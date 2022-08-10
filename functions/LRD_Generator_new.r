library(assertthat)
## Standardized Covariance of Fractional Guassian Process
## r(0) = 1 with r(k) = 0.5*(|k-1|^2H - |k+1|^2H - 2|k|^2H)
FGaussian_cov <- function(n, H){ 
  ####
  # n is the length of ts 
  # H in (0.5, 1) is Hurst index, and alpha = 2 - 2H
  ####
  H2 <- 2*H
  # r:= {r(0), ..., r(n-1)}
  r <- c()
  r[1] <- 2
  k <- 1:(n - 1)
  r[2:n] <- abs(k - 1)^H2 - 2*(k)^H2 + (k + 1)^H2
  # make sure r(0)=1
  r <- r/r[1]
  return(r)
}
# test
if(sys.nframe() == 0){
  FGaussian_cov(10, 0.7)
}

## Standardzied Covariance of FARIMA(0,d,0) Process
## r(0) = 1 with r(k) = 0.5*(|k-1|^2H - |k+1|^2H - 2|k|^2H)
Farima_cov <- function(n, d){
  ####
  # n is the length of ts
  # d in (0, 0.5) is a parameter satisfies alpha = 1 - 2d
  ####
  sigma <- gamma(1- 2*d)/(gamma(1-d)*gamma(d))
  # r:= {r(0), ..., r(n-1)}
  r <- c()
  r[1] <- sigma
  for(k in 1:n){
    r[k+1] <- (k+d-1)/(k-d)*r[k]
  }
  # make sure r(0)=1
  r <- r/r[1]
  return(r[1:n]) 
}

# test
if(sys.nframe() == 0){
  Farima_cov(10, 0.2)
}


#####
# Simulate FGaussian or Farima process via The Davies and Harte method
# by inputting a covariance
LRD_sim <- function(cov){
  n <- length(cov) - 1
  lam <- Re(fft(c(cov[seq(1,n+1, by=1)], rev(cov[c(2:n)])),inverse = T))
  if(!all(lam>=0) ){
    stop('invalid covariance')
  }
  
  W_0N <- rnorm(2)
  V1 <- rnorm(n-1)
  V2 <- rnorm(n-1)
  W <- c(W_0N[1], 1/sqrt(2)*(V1 + V2*1i), W_0N[2], 1/sqrt(2)*rev(V1 - V2*1i))
  ts_sim <- 1/sqrt(2*n)*Re(fft(sqrt(lam)*W, inverse = F))[1:n] #
  return(ts(ts_sim))
}

# Generator of LRD (FGaussian or Farima)
LRD_gen <- function(n, alpha, nrep, option){
  ####
  # n: length of each simulated ts
  # alpha: long-memory exponent
  # nrep: number of generated ts
  # option: 'Farima', 'FGaussian'
  ####
  assert_that(alpha>0 & alpha<1)
  
  if(option == 'FGaussian'){
    H <- 0.5*(2 - alpha)
    cov <- FGaussian_cov(n+1, H)
  }else if(option == 'Farima'){
    d <- 0.5*(1 - alpha)
    cov <- Farima_cov(n+1, d)
  }else{
    stop('Wrong method')
  }
  
  # generate number of nrep ts
  ts_samples <- lapply(c(1:nrep), function(x){
    LRD_sim(cov)
  })
  return(ts_samples)
}

# test
if(sys.nframe() == 0){
  print(LRD_gen(n = 100, alpha = 0.3, nrep = 2, option = 'FGaussian'))
  print(LRD_gen(n = 100, alpha = 0.3, nrep = 2, option = 'Farima'))
  
  # sample autocovariance vs poplutation covariance
  acf(LRD_gen(n = 1000, alpha = 0.3, nrep = 2, option = 'Farima')[[1]])
  lines(c(0:39),Farima_cov(n = 10000, d = 0.35)[1:40])
}



