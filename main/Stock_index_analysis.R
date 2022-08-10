library(tidyverse)
library(EQL)
library(zoo)
library(LongMemoryTS) # local Whittle estimation
library(quantmod) # get data from yahoo finance
library(lubridate) # play with date format
library(aTSA)
library(latex2exp) # latex in ggplot
library(gridExtra) # merging ggplots
library(assertthat)

source("functions/LRD_Generator_new.r")
source("functions/parallel.R")
source("functions/Hermite_rank_estimator.R")

if(sys.nframe() == 0){
  ##### DJI data manipulation
  # get dow jones market data from yahoo finance
  getSymbols("DJI", from="2010-01-01",  src = "yahoo")
  DJI <- DJI %>% 
    dailyReturn() %>% 
    as.data.frame() 
  DJI$Date <- rownames(DJI) %>% as.Date()
  
  # impute missing values for DJI from Marketwatch
  URL <- "https://www.marketwatch.com/investing/index/djia/downloaddatapartial?startdate=05/01/2019%2000:00:00&enddate=05/01/2020%2000:00:00&daterange=d30&frequency=p1d&csvdownload=true&downloadpartial=false&newdates=false"
  DJI_MW <- read.csv(URL, header = T)
  DJI_MW$Date <- as.Date(DJI_MW$Date, format = '%m/%d/%Y')
  DJI_MW <- DJI_MW %>% 
    mutate(across(c(colnames(DJI_MW)[-1]), ~ as.numeric(gsub(',', '', .x))),
           return = (Close - Open)/Open) %>% 
    as_tibble()
  
  # left join two data sets
  DJI <- merge(DJI, DJI_MW, by = 'Date', all = T) %>% 
    dplyr::mutate(daily.returns = ifelse(is.na(daily.returns), return, daily.returns)) %>% 
    dplyr::select(Date, daily.returns) %>% 
    as_tibble() %>% 
    filter(Date < as.Date('2021-01-01')) %>% 
    dplyr::mutate(X_t = abs(daily.returns), Z_t = log(X_t+0.0006))
  
  # test stationarity and normality of Z_t,1
  adf.test(DJI$Z_t-mean(DJI$Z_t), nlag = 50)
  shapiro.test(DJI$Z_t-mean(DJI$Z_t))
  
  #### S&P 500 data manipulation
  SnP500 <- read.csv("data/SnP500_WSJ.csv")
  SnP500 <- SnP500 %>% 
    dplyr::mutate(return = (Close - Open)/Open, Date = as.Date(Date, format = '%m/%d/%y'),
           X_t = abs(return), Z_t =log(X_t+0.0004)) %>% 
    as_tibble()

  ### test stationarity and normality of Z_t,2
  adf.test(SnP500$Z_t-mean(SnP500$Z_t), nlag = 50)
  # shapiro.test(SnP500$Z_t-mean(SnP500$Z_t))
  
  
  ##############
  ############## summary generate function to make input of final plots function
  Hermite_index_analysis_func <- function(df, HE_by, lag=60){
    ####
    # df: a dataframe with Date, X_t, Z_t
    # He_by: by parameter of Hermite rank est for k, DJI=2, Snp500=1
    # lag: the maximum lad for autocorrelation for Z_t
    ####
    
    ### estimate Hermite rank and alpha for the index 
    X_t <- df$X_t
    n <- length(X_t)
    l <- round(2*n^(1/4))
    k_max <- round(3*(n/l)^(1/5))
    k_min <- max(round((k_max)^0.5)+1,2)
    cat('k_min =',k_min, ', k_max =',k_max, '\n')

    est_summary <- HermiteRank_est_func(x = X_t, l = l, 
                                       k_array = seq(k_min, k_max, by = HE_by),
                                       am_option = 'Whittle')
    DJI_alpha_hat <- est_summary$a_hat
    d_hat <- (1-DJI_alpha_hat)/2
    H_hat <- (2-DJI_alpha_hat)/2
    
    ### compute the autocorrelation for the sample and benchmark LRD
    Zt_acf <- acf(df$Z_t, lag.max = lag, type =  "correlation", plot = F)$acf[,1,1]
    cov_comparison_df <- tibble(Lag = seq(0, lag),
                                Emp_Zt = Zt_acf,
                                Farima = FGaussian_cov(n = lag+1, H = H_hat),
                                FGaussian = Farima_cov(n = lag+1, d = d_hat)
    ) %>% 
      gather('type', 'Covariance', -Lag)
    
    return(list(index_df = df, 
                est_summary_df = est_summary,
                cov_comparison_df = cov_comparison_df))
    
  }
  
  
  ######## Final plots function
  stock_plot_func <- function(index_output, name, i, shift){
    ####
    # index_output: output from Hermite_index_analysis_func
    # name <- 'Dow\\,Jones', 'S&P 500'
    # i: subscript, 1 or 2
    # shift: DJI=0.0006, Snp500=0.0004
    ####
    
    ### Visualize dow jones absolute daily return X_t
    Xt_plot <- index_output$index_df %>% 
      ggplot(aes(x = Date, y = X_t)) +
      geom_line() +
      ylab(TeX(paste0('$X_{t,',i,'}$'))) +
      scale_x_continuous(breaks = pretty(DJI$Date, n = 5)) +
      ggtitle(TeX(paste0('$', name,'\\,Index\\,Absolute\\,Daily\\,Return\\,(X_{t,', i,'})$'))) + 
      theme(plot.title = element_text(hjust = 0.5)) 
    
    ### Visualize dow jones absolute daily return X_t
    
    Zt_plot <- index_output$index_df  %>% 
      ggplot(aes(x = Date, y = Z_t - mean(Z_t))) +
      geom_line() +
      ylab(TeX(paste0('$Z_{t,',i,'}$'))) +
      scale_x_continuous(breaks = pretty(DJI$Date, n = 5)) +
      ggtitle(TeX(
        paste0('$Z_{t,',i,'}=log(X_{t,',i,'}+',
               paste0(
                 sprintf( c('%.4f', ')+%.4f$'), 
                          c(shift, abs(round(mean(index_output$index_df$Z_t), digits = 4))))
                 ,collapse = '')
        ))) + 
      theme(plot.title = element_text(hjust = 0.5))
    
    ### histogram for Z_t
    Zt_hist <- ggplot(index_output$index_df, aes(x = Z_t - mean(Z_t))) + 
      geom_histogram(aes(y =..density..),
                     breaks = seq(-3, 3, by = 0.3), 
                     colour = "black", 
                     fill = "grey") +
      stat_function(fun = dnorm, args = list(mean = 0, sd = sd(DJI$Z_t))) +
      xlab(TeX(paste0('$Z_{t,',i,'}'))) +
      ylab('Density') + 
      ggtitle(TeX(paste0('$Histogram\\,of\\,Z_{t,',i,'$'))) + 
      theme(plot.title = element_text(hjust = 0.5))
    
    ### Comparison between sample covariance and LRD covariance
    cov_plot <- index_output$cov_comparison_df %>% 
      ggplot(aes(x = Lag, y = Covariance, linetype = as.factor(type)) ) +
      geom_line() +
      scale_linetype_discrete(labels=lapply(c('$Sample$', 'Farima', 'FGaussian'), TeX)) +
      ggtitle(TeX(
        paste0('$Sample\\,Correlation\\,for\\,Z_{t,',i,'}\\,vs\\,LRD\\,Correlation\\,(\\alpha =$', 
               round(index_output$est_summary_df[4], digits = 4), ')')
        )) + 
      ylab(label = 'Correlation')+
      guides(linetype=guide_legend(title=NULL)) +
      theme( legend.background = element_blank(),
             legend.key = element_blank(),
             legend.position= c(0.85, 0.8),
             legend.box.background = element_rect(colour = "black"),
             plot.title = element_text(hjust = 0.5))
    
    Merged_plots <- grid.arrange(Xt_plot, Zt_plot, Zt_hist, cov_plot, nrow = 2)
    
  }
  
  #### making final plot for DJI and S&p 500
  DJI_index_output <- Hermite_index_analysis_func(df = DJI, 
                                              HE_by=2, lag=60)
  DJI_plots <- stock_plot_func(index_output = DJI_index_output,
                  name = 'Dow\\,Jones',
                  i=1,
                  shift = 0.0006)
  
  SnP500_index_output <- Hermite_index_analysis_func(df = SnP500, 
                                                  HE_by=1, lag=60)
  SnP500_plots <- stock_plot_func(index_output = SnP500_index_output,
                               name = 'S&P\\,500',
                               i=2,
                               shift = 0.0004)
  # ggsave("DJI_plots.png", plot = DJI_plots,  width = 11, height = 7)
  # ggsave("SnP500_plots.png", plot = SnP500_plots,  width = 11, height = 7)
}





# ##### Hermite rank estimation of absolute daily return for Index data
# 
# Estimate_m_emp_func <- function(x, l, k_array, am_option){
#   ####
#   # x: a single time series
#   # l: block size
#   # k_array: different k for squared blocks
#   # am_option: 'Whittle', 'periodogram'
#   ####
#   
#   assert_that(am_option %in% c('Whittle', 'periodogram'))
#   
#   ### estimate 2a by covariance method
#   # \hat{alpha} = [log(cov(k1)) - log(cov(k2))]/[2(log(k2) - log(k1))]
#   # x: a single ts with length n
#   # estimate alpha*m 
#   n <- length(x)
#   am_est_func <- function(x, n){
#     if(am_option == 'Whittle'){
#       am_hat <- (1 - 2 *GSE(x,  m=floor(1+n^0.65)))
#     }else if(am_option == 'periodogram'){
#       am_hat <- (1 - 2 *gph(x,  m=floor(1+n^0.7)))
#     }
#     return(am_hat)
#   }
#   am_hat <- am_est_func(x=x, n=n) %>%
#     abs() %>% 
#     min(0.99) %>% 
#     max(0.01)
#   
#   # V_hat = V_hat/\ell^{\am}
#   x_bar_mv <- rollmean(x, l)
#   x_bar_2 <- (x_bar_mv - mean(x_bar_mv))^2
#   V_hat <- mean(x_bar_2)
#   
#   # log_cov_array return length(k_array) of log(cov_k)
#   Nl <- length(x_bar_2)
#   if(k_array[length(k_array)]*l >= Nl){
#     stop('kl > Nl')
#   }
#   log_cov_array <- abs(sapply(k_array, function(k){
#     mean(x_bar_2[seq(1, Nl-k*l)]*x_bar_2[seq(1+k*l, Nl)]) - V_hat^2
#   })) %>% 
#     log()
#   
#   log_k_array <- log(k_array)
#   a_hat <- -0.5*lm(log_cov_array~log_k_array, 
#                    data = data.frame(log_cov_array, log_k_array))$coefficients[2]
#   a_hat <- a_hat %>% abs()
# 
#   m_hat <- am_hat/a_hat 
#   
#   return(data.frame(m_hat = m_hat, am_hat = am_hat, a_hat = a_hat))
# }


