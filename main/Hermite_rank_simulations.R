library(tidyverse)
library(EQL)
library(zoo)
library(LongMemoryTS) # local Whittle estimation
library(quantmod) # get data from yahoo finance
library(lubridate) # play with date format
library(aTSA)
library(latex2exp) # latex in ggplot
library(gridExtra) # merging ggplots
library(assertthat) # check function variable conditions
library(pbapply) # progress bar
library(knitr) # generate latex code for df

source("functions/LRD_Generator_new.r")
source("functions/parallel.R")


#### function for estimate hermite rank by covariance method in smulation study
Estimate_m_cov_Sim <- function(n, l, alpha, k_array, m_array, 
                               nrep, G_list, option, am_option){
  ####
  # n: length of each ts
  # l: block size
  # alpha: long memory exponent
  # k_array: different k for squared blocks
  # m_array: the hermite rank of G_list
  # nrep: number of LRD series to generate
  # G_list: a list of functions G(.)
  # gen_option: 'Farima' or 'FGaussian'
  # am_option: 'Whittle' or 'periodogram'
  ####
  
  assert_that(option %in% c('Farima', 'FGaussian'))
  assert_that(am_option %in% c('Whittle', 'periodogram'))
  assert_that(length(G_list) == length(m_array))
  assert_that(k_array[length(k_array)]*l <= (n-l+1))
  
  # estimate function for alpha*m
  am_est_func <- function(x, n){
    if(am_option == 'Whittle'){
      am_hat <- (1 - 2 *GSE(x,  m=floor(1+n^0.65)))
    }else if(am_option == 'periodogram'){
      am_hat <- (1 - 2 *gph(x,  m=floor(1+n^0.7)))
    }
    return(am_hat)
  }
  
  ### estimate 2a by covariance method
  # \hat{alpha} = [log(cov(k1)) - log(cov(k2))]/[2(log(k2) - log(k1))]
  # m_hat_mat: length(nrep) by length(m_array) matrix
  summary_table_nrep <- lapply(c(1:nrep), function(i){ 
    repeat{
      # x: a single ts with length n
      zt <- LRD_gen(n = n, 
                    alpha = alpha, 
                    nrep = 1, 
                    option = option)[[1]]
      
      Gx_list <- lapply(seq(1, length(G_list)), function(j){
        G_list[[j]](zt)
      })
      
      am_hat_list <- sapply(Gx_list, function(Gx){
        am_hat <- am_est_func(x = Gx, n = n)
        keep_bool <- (am_hat>0 & am_hat< 1)
        return(c(am_hat, keep_bool))
      }) %>% 
        t() %>% 
        as_tibble() %>% 
        `colnames<-`(c('am_hat', 'keep_bool'))
      
      # check the generated ts quality
      if(sum(am_hat_list$keep_bool) == length(G_list)){
        break
      }
    }

    summary_table <- sapply(seq(1, length(G_list)), function(j){
      m <- m_array[j]
      x <- Gx_list[[j]]
      am_hat <- am_hat_list$am_hat[j]
      
      # V_hat = V_hat/\ell^{\am}
      x_bar_mv <- rollmean(x, l)
      x_bar_2 <- (x_bar_mv - mean(x_bar_mv))^2
      V_hat <- mean(x_bar_2)
      
      # log_cov_array: vector of length(k_array) of log(cov_k)
      Nl <- length(x_bar_2)
      log_cov_array <- abs(sapply(k_array, function(k){
        mean(x_bar_2[seq(1, Nl-k*l)]*x_bar_2[seq(1+k*l, Nl)]) - V_hat^2
      })) %>% 
        log()
      
      log_k_array <- log(k_array)
      a_hat <- -0.5*lm(log_cov_array~log_k_array, 
                       data = data.frame(log_cov_array, log_k_array))$coefficients[2]
      
      a_hat <- a_hat %>% 
        abs() %>% 
        min(1) %>% 
        max(0)
      
      m_hat <- max(round(am_hat/a_hat), 1) 
      
      a_adj_hat <- am_hat/m_hat
      return(c(m, m_hat, a_hat, a_adj_hat, am_hat))
    }) %>%
      t() %>%
      as_tibble() %>%
      `colnames<-`(c('m', 'm_hat', 'a_hat', 'a_adj_hat', 'am_hat'))
    return(summary_table)
  }) %>% 
    Reduce(f = rbind)
  
  return(summary_table_nrep)
}


#########################
#########################


### m_hat estimation summary for Whittle and Log-periodogram

if(sys.nframe() == 0){
  #### Specify a list of transformation function G(.)
  G_list <- list()
  G_list[[1]] <- function(x) {
    sin(x)
  }
  G_list[[2]] <- function(x) {
    hermite(x, n = 2, prob = T) + 
      hermite(x, n = 4, prob = T) +
      hermite(x, n = 5, prob = T)
  }
  G_list[[3]] <- function(x) {
    hermite(x, n = 3, prob = T) + hermite(x, n = 4, prob = T)
  }
  
  n0_array <- c(800, 2000, 5000)
  nrep <- 1000
  set.seed(44)
  
  Hermite_sim_result_ss <- mclapply.hack(c('Whittle', 'periodogram'), function(am_method){
    ### output for alpha=0.3
    Hermite_sim_a_0.3 <- lapply(n0_array, function(n0){
      l0 <- 2*n0^(1/4)
      k_max <- 3*(n0/l0)^(1/5)
      k_range <- seq(round(max((k_max)^0.5,2)),round(k_max))
      base::print(paste('n0=', n0, 'l0=', l0, 'k_max=', k_max, 'k_array=', c(round(max((k_max)^0.5,2)),round(k_max))))
      m_est_result <- Estimate_m_cov_Sim(n = n0, l = round(l0), alpha = 0.3, nrep = nrep,
                                         k_array = k_range, m_array = c(1, 2, 3),
                                         G_list = G_list, option = 'FGaussian',
                                         am_option = am_method)
      m_hat_summary <- m_est_result %>% 
        group_split(m) %>% 
        lapply(function(df){
          table(df$m_hat)
        })
      
      a_hat_summary <- m_est_result %>% 
        group_split(m) %>% 
        lapply(function(df){
          df$a_hat
        })
      
      a_adj_hat_summary <- m_est_result %>% 
        group_split(m) %>% 
        lapply(function(df){
          df$a_adj_hat
        })
      
      return(list(m_hat_summary = m_hat_summary,
                  a_hat_summary = a_hat_summary,
                  a_adj_hat_summary = a_adj_hat_summary))
    })
    
    ### output for alpha=0.15
    Hermite_sim_a_0.15 <- lapply(n0_array, function(n0){
      l0 <- 2*n0^(1/4)
      k_max <- 3*(n0/l0)^(1/5)
      k_range <- seq(round(max((k_max)^0.5,2)),round(k_max))
      base::print(paste('n0=', n0, 'l0=', l0, 'k_max=', k_max, 'k_array=', c(round(max((k_max)^0.5,2)),round(k_max))))
      m_est_result <- Estimate_m_cov_Sim(n = n0, l = round(l0), alpha = 0.15, nrep = nrep,
                                         k_array = k_range, m_array = c(1, 2, 3),
                                         G_list = G_list, option = 'FGaussian',
                                         am_option = am_method)
      m_hat_summary <- m_est_result %>% 
        group_split(m) %>% 
        lapply(function(df){
          table(df$m_hat)
        })
      
      a_hat_summary <- m_est_result %>% 
        group_split(m) %>% 
        lapply(function(df){
          df$a_hat
        })
      
      a_adj_hat_summary <- m_est_result %>% 
        group_split(m) %>% 
        lapply(function(df){
          df$a_adj_hat
        })
      
      return(list(m_hat_summary = m_hat_summary,
                  a_hat_summary = a_hat_summary,
                  a_adj_hat_summary = a_adj_hat_summary))
    })
    
    return(list(Hermite_sim_a_0.3 = Hermite_sim_a_0.3,
                Hermite_sim_a_0.15 = Hermite_sim_a_0.15))
  })
  
  names(Hermite_sim_result_ss) <- c('Whittle', 'periodogram')
  
  #saveRDS(Hermite_sim_result_ss, 'Hermite_sim_result_ss.rds')
  
  if(exists('Hermite_sim_result_ss')){
    Hermite_sim_result <- Hermite_sim_result_ss
  }else{
    Hermite_sim_result <- readRDS("C:\\Users\\howard\\Desktop\\PhD Research\\LRD Resampling\\Proj 2 Estimation of Hermite rank\\Hermite_rank_code\\cache\\Hermite_sim_result_ss.rds")
  }
  
  ### remark: whittle result based on Hermite_sim_result
  ###         periodogram result based on Hermite_sim_result
  ### latex code for estimation of m
  
  mhat_latex <- lapply(Hermite_sim_result, function(lst_bymethod){
    lapply(lst_bymethod, function(lst_by_alpha){
      lapply(c(1:3), function(i){
        # each i for n = 800, 2000, 5000
        lst <- lst_by_alpha[[i]]$m_hat_summary
        m_table_latex <- sapply(lst, function(lst_m){
          lst <- c(lst_m[1:4], sum(lst_m[5:length(lst)]))
        }) %>% t() %>% 
          as_tibble() %>% 
          knitr::kable(format = 'latex')
      }) %>% 
        `names<-`(paste0('n=', c(800,2000,5000)))
    })
  })
  
  mhat_latex$periodogram$Hermite_sim_a_0.3$`n=800`
  mhat_latex$periodogram$Hermite_sim_a_0.3$`n=2000`
  mhat_latex$periodogram$Hermite_sim_a_0.3$`n=5000`
  
  
  lapply(c(1:3), function(i){
    # each i for n = 800, 2000, 5000
    lst <- Hermite_sim_result$periodogram$Hermite_sim_a_0.3[[i]]$m_hat_summary
    m_table_latex <- sapply(lst, function(lst_m){
      lst <- c(lst_m[1:4], sum(lst_m[5:length(lst)]))
    }) %>% t() %>% 
      as_tibble() %>% 
      knitr::kable(format = 'latex')
  }) %>% 
    `names<-`(paste0('n=', c(800,2000,5000)))
  
}

