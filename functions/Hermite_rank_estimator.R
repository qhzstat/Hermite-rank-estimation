library(tidyverse)
library(LongMemoryTS) # for am_hat

HermiteRank_est_func <- function(x, l, k_array, am_option){
  ####
  # x: a single time series
  # l: block size
  # k_array: different lag k for squared blocks
  # am_option: 'Whittle', 'periodogram'
  ####
  
  assert_that(length(k_array)>=2)
  assert_that(am_option %in% c('Whittle', 'periodogram'))

  # estimate alpha*m
  n <- length(x)
  am_est_func <- function(x, n){
    if(am_option == 'Whittle'){
      am_hat <- (1 - 2 *GSE(x,  m=floor(1+n^0.65)))
    }else if(am_option == 'periodogram'){
      am_hat <- (1 - 2 *gph(x,  m=floor(1+n^0.7)))
    }
    return(am_hat)
  }
  am_hat <- am_est_func(x=x, n=n)
  
  if(am_hat<0 | am_hat >1){
    warning('am_hat estimation out of (0, 1)')
  }
  am_hat <- abs(am_hat)
  
  
  ### estimate 2a by covariance method
  # \hat{alpha} = [log(cov(k1)) - log(cov(k2))]/[2(log(k2) - log(k1))]
  x_bar_mv <- rollmean(x, l)
  x_bar_2 <- (x_bar_mv - mean(x_bar_mv))^2
  V_hat <- mean(x_bar_2)
  
  # log_cov_array return length(k_array) of log(cov_k)
  Nl <- length(x_bar_2)
  if(k_array[length(k_array)]*l >= Nl){
    stop('kl > Nl')
  }
  log_cov_array <- abs(sapply(k_array, function(k){
    mean(x_bar_2[seq(1, Nl-k*l)]*x_bar_2[seq(1+k*l, Nl)]) - V_hat^2
  })) %>% 
    log()
  
  log_k_array <- log(k_array)
  a_hat <- -0.5*lm(log_cov_array~log_k_array, 
                   data = data.frame(log_cov_array, log_k_array))$coefficients[2]
  a_hat <- a_hat %>% abs()
  
  m_hat <- am_hat/a_hat 
  
  m_int_hat <- max(1, round(m_hat))
  
  return(data.frame(m_hat = m_hat, m_int_hat = m_int_hat,
                    am_hat = am_hat, a_hat = a_hat))
}


### test
if(sys.nframe() == 0){
  ts_list <- LRD_gen(n = 1000, alpha = 0.1, nrep = 5, option = 'FGaussian')
  
  lapply(ts_list, function(ts){
    n <- length(ts) 
    l <- round(2*n^(1/4))
    k_max <- 3*(n/l)^(1/5)
    k_min <- max(2, round(k_max^(1/2)))
    k_array <- seq(k_min, max(k_min+1, round(k_max)))
    HermiteRank_est_func(x = ts, l = l, k_array = k_array, am_option = 'Whittle')
  })
  
  lapply(ts_list, function(ts){
    n <- length(ts) 
    l <- round(2*n^(1/4))
    k_max <- 3*(n/l)^(1/5)
    k_min <- max(2, round(k_max^(1/2)))
    k_array <- seq(k_min, max(k_min+1, round(k_max)))
    HermiteRank_est_func(x = ts, l = l, k_array = k_array, am_option = 'periodogram')
  })
    
}








