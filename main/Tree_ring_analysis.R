library(dplR) # read for tree ring data
library(tidyverse)
library(zoo)
library(LongMemoryTS) # local Whittle estimation
library(pbapply) # progress bar
library(latex2exp) # latex in ggplot
library(gridExtra) # merging ggplots

source("functions/LRD_Generator_new.r")
source("functions/parallel.R")
source("functions/Hermite_rank_estimator.R")

path <- getwd()
path_crn <- paste0(path, '/data/crn_data')

### read all readable crn files with length(ts)>500
file_names <- list.files(path = path_crn, pattern = '.crn')
crn_list <- lapply(file_names, function(name){
  tryCatch(({
    df <- dplR::read.crn(fname = name, header = T) %>% 
      `colnames<-`(c('treering', 'samp.depth'))
  }), error = function(e){
    cat('wrong crn file')
    df <- NULL
  }, finally = ({
    if(is.data.frame(df)){
      if(dim(df)[1] <= 500){
        df <- NULL
      }
    }
    df
  }))
})
crn_list <- crn_list[-which(sapply(crn_list, is.null))]
#saveRDS(object = crn_list, file = "cache/crn_list.rds")

##### function for tree rings filter
keep_ts <- function(x){
  x <- na.omit(x)
  n <- length(x)
  am_hat <- (1 - 2 *GSE(x,  m=floor(1+n^0.65)))
  ### check am_hat in (0, 1), o.w drop the sample
  keep <- (am_hat>= 0.05 & am_hat <= 0.95)
  return(c(keep, am_hat))
}


##### Hermite rank estimation of tree ring ts

Estimate_m_tree_func <- function(x, l, k_array, am_option){
  ####
  # x: a single time series
  # l: block size
  # k_array: different k for squared blocks
  ####
  
  assert_that(am_option %in% c('Whittle', 'periodogram'))
  assert_that(length(k_array) >= 2)
  
  ### estimate alpha*m via Local Whittle estimator
  x <- na.omit(x)
  n <- length(x)
  if(am_option == 'Whittle'){
    am_hat <- (1 - 2 *GSE(x,  m=floor(1+n^0.65)))
  }else if(am_option == 'periodogram'){
    am_hat <- (1 - 2 *gph(x,  m=floor(1+n^0.7)))
  }else{
    stop('invalid am estimation method')
  }
  
  ### check am_hat in (0, 1), o.w drop the sample
  keep <- (am_hat>= 0.05 & am_hat <= 0.95)
  if(keep == F){
    return(list(keep = keep, est = NULL))
  }
  
  estimation_df <- HermiteRank_est_func(x = x, l=l, k_array = k_array, am_option = am_option)
  
  return(list(keep = keep, 
              est = data.frame(n=n, estimation_df)))
}


### function for generate FGaussian ts and estimate Hermite rank
tree_sim <- function(n, alpha, nrep, option, am_option){
  ####
  # n: length of ts
  # alpha: long-momory exponent
  # nrep: number of replications to generate
  # am_option: 'Whittle', 'periodogram'
  ####
  
  assert_that(am_option %in% c('Whittle', 'periodogram'))
  
  x_list <- LRD_gen(n=n, alpha = alpha, nrep = nrep, option = option)
  x2_list <- lapply(x_list, function(x) x^2)
  
  # parameters set up for Hermite rank estimation
  l <- round(2*n^(1/4))
  k_max <- 3*(n/l)^(1/5)
  k_min <- round(max((k_max)^0.5,2))
  k_max <- max(k_min+1, round(k_max))
  k_array <- seq(k_min, k_max, by = 1)
  
  # get Hermite rank estimation for nrep of x and x^2
  x_m_hat <- sapply(x_list, function(x){
    m_est <- HermiteRank_est_func(x=x, l=l, k_array = k_array, am_option = am_option)$m_int_hat
  })
  x2_m_hat <- sapply(x2_list, function(x){
    m_est <- HermiteRank_est_func(x=x, l=l, k_array = k_array, am_option = am_option)$m_int_hat
  })
  
  return(tibble(x_m_hat=x_m_hat, x2_m_hat=x2_m_hat))
  
}

if(sys.nframe() == 0){
  # read cleaned crn data from cache if crn_list does not exist
  if(exists('crn_list')== F){
    crn_list <- readRDS("cache\\crn_list.rds")
  }
  
  # histogram of length n for all filtered ts by local Whittle estimation
  tree_m_stats <- mclapply(crn_list, function(crn){
    Xt <- na.omit(crn$treering - mean(crn$treering))
    if(length(Xt) <= 500){
      return(NULL)
    }
    keep_Xt <- keep_ts(Xt)
    if(keep_Xt[1] & keep_ts(Xt^2)[1]){
      n <- length(Xt)
      am_hat <- keep_Xt[2]
      return(c(n, am_hat))
    }else{
      return(NULL)
    }
  })  %>% 
    Reduce(f=rbind) %>% 
    as_tibble() %>% 
    `colnames<-`(c('Length', 'am_hat'))
  
  # 5 number summary 
  tree_m_stats$Length %>% summary()
  tree_m_stats$am_hat %>% summary()
  
  # histogram for length
  tr_length_hist <- tree_m_stats %>%
    ggplot(aes(x=Length)) +
    geom_histogram(aes(y=..density..), color='black', fill=c("#69b3a2"), position = 'identity', binwidth = 100) +
    xlim(c(0, 8000)) +
    xlab('Length') +
    ylab('Frequency') + 
    ggtitle('Histogram of Sample Sizes') + 
    theme(plot.title = element_text(hjust = 0.5))
  
  # boxplot for length
  tr_length_boxplot <- ggplot(data=tree_m_stats[tree_m_stats$Length<=1500, ], aes(x='', y=Length)) +
    geom_boxplot(color='black', fill=c("#69b3a2")) +
    xlab('') +
    ylab('Length') +
    ggtitle('Boxplot of Sample Sizes') + 
    theme(plot.title = element_text(hjust = 0.5))
  
  # histogram for am_hat
  tr_am_hist <- tree_m_stats %>%
    ggplot(aes(x=am_hat)) +
    geom_histogram(aes(y=..density..), color='black', fill=c("#69b3a2"), position = 'identity') +
    xlab(TeX('$\\widehat{\\alpha m}$')) + 
    ylab('Frequency') +
    ggtitle('Histogram of Memory Parameter Estimates') + 
    theme(plot.title = element_text(hjust = 0.5))
  
  # boxplot for am_hat
  tr_am_boxplot <- ggplot(data=tree_m_stats, aes(x='', y=am_hat)) +
    geom_boxplot(color='black', fill=c("#69b3a2")) +
    xlab('') +
    ylab(TeX('$\\widehat{\\alpha m}$')) +
    ylim(c(0,1)) +
    ggtitle('Boxplot of Memory Parameter Estimates') + 
    theme(plot.title = element_text(hjust = 0.5)) 
  
  
  set.seed(44)
  nrep <- 100
  
  
  tree_est_list <- mclapply.hack(c('Whittle', 'periodogram'), function(method){
    tree_est_bymethod <- lapply(crn_list, function(crn){
      Xt <- na.omit(crn$treering - mean(crn$treering))
      if(length(Xt) <= 500){
        return(NULL)
      }
      
      # set up parameters form hermite rank estimation
      n0 <- length(Xt)
      l0 <- round(2*n0^(1/4))
      k_max <- 3*(n0/l0)^(1/5)
      k_min <- round(max((k_max)^0.5,2))
      k_max <- max(k_min+1, round(k_max))
      k_array <- seq(k_min, k_max, by = 1)
      Xt_est <- Estimate_m_tree_func(x=Xt, l =l0, k_array = k_array, 
                                     am_option = method)
      if(Xt_est$keep == F){
        return(NULL)
      }else{
        Xt2_est <- Estimate_m_tree_func(x=Xt^2, l =l0, k_array = k_array, 
                                        am_option = method)
        if(Xt2_est$keep == F){
          return(NULL)
        }
        # estimation summary table for Xt and Xt^2
        est_df <- bind_rows(Xt_est$est , Xt2_est$est) %>% 
          as_tibble() %>% 
          dplyr::mutate(Xt=c(1,2))
      }
      return(est_df)
    })
    
    # remove null element 
    tree_est_bymethod <- tree_est_bymethod[-which(sapply(tree_est_bymethod,is.null))]
    
    tbl <- tree_est_bymethod[[1]]
    # matrix for m_hat of Xt and Xt^2
    tree_output_list <- tree_est_bymethod %>% lapply(function(tbl){
      # table of m_int_hat of tree ring data
      m_int_hat <- tbl$m_int_hat # vector with length 2
      
      # table of m_hat of generated Gaussian Xt and Xt^2
      par_set <- tbl %>% 
        filter(Xt==1)
      
      sim_est <- tree_sim(n= par_set$n, alpha = par_set$am_hat, nrep = nrep, 
                          option = 'FGaussian', am_option = method)
      
      return(list(m_int_hat = m_int_hat, sim_est = sim_est))
    }) 
    
  })
  names(tree_est_list) <- c('Whittle', 'periodogram') 
  #saveRDS(object = tree_est_list, file = "C:\\Users\\howard\\Desktop\\PhD Research\\LRD Resampling\\Proj 2 Estimation of Hermite rank\\Code\\cache\\tree_est_list.rds")
  
  
  tree_est_list <- readRDS(file = "cache\\tree_est_list.rds")
  ### Hermite Hermite Rank Estimation Comparison Plot for Xt vs Xt^2, null for Whittle
  plot_names <- c('Samplewise Comparison from Tree Rings',
                  'Samplewise Comparison (Log-periodogram)')
  hist_names <- c('Histogram Comparison from Tree Rings',
                  'Histogram Comparison from Tree Rings (Log-periodogram)')
  sim_names <- c('Histogram Comparison from Simulation',
                 'Histogram Comparison from Simulation (Log-periodogram)')
  
  ### for two different am estimation methods: Whittle, log-periodogram
  Treering_plot_list <- list()
  tree_sim_comparison <- list()
  for(i in c(1,2)){
    tree_sim_comparison[[i]] <- list()
    ### count the frequency of m_hat(Xt) > m_hat(Xt^2)
    tree_sim_comparison[[i]]$tree_m_comparison <- sapply(tree_est_list[[i]], function(lst){
      lst$m_int_hat
    }) %>% 
      t() %>% 
      as_tibble() %>% 
      `colnames<-`(c('Xt', 'Xt2')) %>% 
      apply(1, function(x){
        if(x[2]>x[1]){
          return(1)
        }else if(x[2] == x[1]){
          return(0)
        }else{
          return(-1)
        }
      }) %>% 
      table()
    
    # m_comparison for simulations
    tree_sim_comparison[[i]]$sim_m_comparison <- sapply(tree_est_list[[i]], function(lst){
      lst$sim_est %>% 
        apply(1, function(x){
          if(x[2]>x[1]){
            return(1)
          }else if(x[2] == x[1]){
            return(0)
          }else{
            return(-1)
          }
        }) %>% 
        table()
    }) %>% 
      t() %>% 
      colSums()
    
    ### plots for Hermite rank estimation of tree rings data
    # remark: type 1 for Xt and type 2 for Xt^2
    tree_m_hat_tbl <- sapply(tree_est_list[[i]], function(lst){
      lst$m_int_hat
    }) %>% 
      t() %>% 
      as_tibble() %>% 
      `colnames<-`(c(1,2)) %>% 
      # arrange(`1`) %>% #
      pivot_longer(cols = everything(), names_to = 'type', values_to = 'm_hat') %>% 
      dplyr::mutate(m_hat = ifelse(m_hat>20, 20, m_hat), 
                    index = rep.int(rep(seq(1, 0.5*n())), 
                                    times = rep(2, times= 0.5*n())))
    
    Treering_plot_list[[i]] <- list()
    Treering_plot_list[[i]][[1]] <- ggplot(tree_m_hat_tbl, aes(x = index, y = m_hat, colour = type)) +
      geom_point() + 
      scale_colour_manual(name = '',#TeX('$G(X_t)$'),
                          values=c("#69b3a2", "#404080"), 
                          labels= lapply(c('$X_t$', '$X_t^2$'), TeX)) + 
      scale_y_continuous(breaks = seq(0, 20, by = 5), lim = c(0, 20), label = c(seq(0, 15, by = 5), '20+')) +
      ylab(TeX('$\\widehat{m}$')) +
      ggtitle(plot_names[i]) + 
      theme(plot.title = element_text(hjust = 0.5)) 
    
    tree_m_hat_tbl <- tree_m_hat_tbl %>% 
      dplyr::mutate(m_hat_trunc = ifelse(m_hat>=10, 10, m_hat))
    
    Treering_plot_list[[i]][[2]] <- ggplot(tree_m_hat_tbl , aes(x=m_hat_trunc, fill=type)) +
      geom_histogram(aes(y=..density..),  alpha=0.4, position = 'identity', binwidth = 1) +
      guides(fill=guide_legend(title='')) + #TeX('$G(X_t)$') 
      scale_fill_manual(values=c("#69b3a2", "#404080"), labels=lapply(c('$X_t$', '$X_t^2$'), TeX)) +
      scale_x_continuous(breaks = seq(1, 10), lim = c(0, 11), labels = c(1:9, '9+')) +
      xlab(TeX('$\\widehat{m}$')) +
      ylab('Density') + 
      ylim(c(0, 0.6)) + 
      ggtitle(hist_names[i]) + 
      theme(plot.title = element_text(hjust = 0.5))
    
    ### plots for Hermite rank estimation of simulated FGuassian
    sim_m_hat_tbl <- lapply(tree_est_list[[i]], function(lst){
      lst$sim_est
    }) %>% 
      Reduce(f = rbind) %>% 
      apply(2, function(m){
        ifelse(m>=10, 10, m)
      }) %>% 
      as_tibble() %>% 
      `colnames<-`(c(1,2)) %>% 
      pivot_longer(cols = everything(), names_to = 'type', values_to = 'm_hat')
    
    Treering_plot_list[[i]][[3]] <- ggplot(sim_m_hat_tbl , aes(x=m_hat, fill=type)) +
      geom_histogram(aes(y=..density..),  alpha=0.4, position = 'identity', binwidth = 1) +
      guides(fill=guide_legend(title='')) + #TeX('$G(X_t)$')
      scale_fill_manual(values=c("#69b3a2", "#404080"), labels=lapply(c('$Z_t^*$', '$(Z_t^*)^2$'), TeX)) +
      scale_x_continuous(breaks = seq(1, 10), lim = c(0, 11), labels = c(1:9, '9+')) +
      ylab('Density') + 
      ylim(c(0, 0.6)) + 
      # scale_y_continuous(breaks = c(0, 10000, 20000, 30000), lim = c(0, 32000)) +
      xlab(TeX('$\\widehat{m}$')) +
      ggtitle(sim_names[i]) + 
      theme(plot.title = element_text(hjust = 0.5))
    
    ### Chi-square test for comparing two histograms
    sim_m_hat_freq <- sim_m_hat_tbl %>% 
      group_by(type, m_hat) %>% 
      summarise(n = n()) %>% 
      mutate(freq = n/sum(n))
    
    tree_m_hat_freq <- tree_m_hat_tbl %>% 
      group_by(type, m_hat_trunc) %>% 
      summarise(n = n()) %>% 
      mutate(freq = n/sum(n))
    
    tree_sim_comparison[[i]]$Chisq_pvalues <- sapply(c(1, 2), function(j){
      sim_subset <- sim_m_hat_freq %>% 
        filter(type == j) %>% 
        select(freq, n)
      tree_subset <- tree_m_hat_freq %>% 
        filter(type == j) %>% 
        select(freq, n)
      
      # chi-square test between hist of Simulations vs tree ring for Xt and Xt^2
      chi_stat <- sum((sim_subset$freq - tree_subset$freq)^2/
                        (sim_subset$freq/sum(sim_subset$n) + tree_subset$freq/sum(tree_subset$n)))
      p_value <- 1- pchisq(chi_stat, df = 9, lower.tail = T)
    }) %>% 
      `names<-`(c('Xt','Xt2'))
  }
  
  # Bubble plot for the proportion of tree ring data for (m^(1), m^(2)) 
  tree_ranks_list <- tree_est_list$Whittle
  proportion_plot <- lapply(tree_ranks_list, function(lst){
    # select the results by Whittle estimator
    lst[[1]] 
  }) %>%
    Reduce(f=rbind) %>%
    pmin(5) %>%
    as_tibble() %>%
    `colnames<-`(c('x1', 'x2')) %>%
    group_by(x1, x2) %>%
    summarise(n=n()) %>%
    ungroup() %>%
    mutate(proportion = n/sum(n)) %>%
    ggplot(aes(x=factor(x1), y=factor(x2), size = proportion)) +
    geom_point(alpha=0.7) +
    scale_x_discrete(
      labels = c('1'='1', '2'='2', '3'='3', '4'='4','5'='4+'
      )
    ) +
    scale_y_discrete(
      labels = c('1'='1', '2'='2', '3'='3', '4'='4','5'='4+'
      )
    ) +
    scale_size(range = c(3, 10), name="Proportion") +
    xlab(TeX('$\\widehat{m}^{(1)}$')) +
    ylab(TeX('$\\widehat{m}^{(2)}$')) +
    ggtitle(TeX('Pairwise Hermite Rank from Tree Rings')) +
    theme(plot.title = element_text(hjust = 0.5))
  
  tr_sample_plots <- grid.arrange(Treering_plot_list[[1]][[1]], proportion_plot, nrow = 1)
  
  tr_hists <- grid.arrange(Treering_plot_list[[1]][[2]], Treering_plot_list[[1]][[3]], nrow = 1)
  
  tr_length <- grid.arrange(tr_length_boxplot, tr_length_hist, nrow = 1)
  
  tr_am <- grid.arrange(tr_am_boxplot, tr_am_hist, nrow = 1)
  
  ggsave("plots/TreeRing_length_plot.png", plot = tr_length,  width = 11, height = 3.5)
  ggsave("plots/TreeRing_am_plot.png", plot = tr_am,  width = 11, height = 3.5)
  ggsave("plots/TreeRing_scatter_plot.png", plot = tr_sample_plots,  width = 11, height = 3.5)
  ggsave("plots/TreeRing_hist_plot.png", plot = tr_hists,  width = 11, height = 3.5)
}


# if(sys.nframe() == 0){
#   # read cleaned crn data from cache if crn_list does not exist
#   if(exists('crn_list')== F){
#     crn_list <- readRDS("cache\\crn_list.rds")
#   }
#   
#   
#   set.seed(44)
#   nrep <- 100
#   
#   tree_est_list <- mclapply.hack(c('Whittle', 'periodogram'), function(method){
#     tree_est_bymethod <- lapply(crn_list, function(crn){
#       Xt <- na.omit(crn$treering - mean(crn$treering))
#       if(length(Xt) <= 500){
#         return(NULL)
#       }
#       
#       # set up parameters form hermite rank estimation
#       n0 <- length(Xt)
#       l0 <- round(2*n0^(1/4))
#       k_max <- 3*(n0/l0)^(1/5)
#       k_min <- round(max((k_max)^0.5,2))
#       k_max <- max(k_min+1, round(k_max))
#       k_array <- seq(k_min, k_max, by = 1)
#       Xt_est <- Estimate_m_tree_func(x=Xt, l =l0, k_array = k_array, 
#                                      am_option = method)
#       if(Xt_est$keep == F){
#         return(NULL)
#       }else{
#         Xt2_est <- Estimate_m_tree_func(x=Xt^2, l =l0, k_array = k_array, 
#                                         am_option = method)
#         if(Xt2_est$keep == F){
#           return(NULL)
#         }
#         # estimation summary table for Xt and Xt^2
#         est_df <- bind_rows(Xt_est$est , Xt2_est$est) %>% 
#           as_tibble() %>% 
#           dplyr::mutate(Xt=c(1,2))
#       }
#       return(est_df)
#     })
#     
#     # remove null element 
#     tree_est_bymethod <- tree_est_bymethod[-which(sapply(tree_est_bymethod,is.null))]
#     
#     tbl <- tree_est_bymethod[[1]]
#     # matrix for m_hat of Xt and Xt^2
#     tree_output_list <- tree_est_bymethod %>% lapply(function(tbl){
#       # table of m_int_hat of tree ring data
#       m_int_hat <- tbl$m_int_hat # vector with length 2
#       
#       # table of m_hat of generated Gaussian Xt and Xt^2
#       par_set <- tbl %>% 
#         filter(Xt==1)
#       
#       sim_est <- tree_sim(n= par_set$n, alpha = par_set$am_hat, nrep = nrep, 
#                option = 'FGaussian', am_option = method)
#       
#       return(list(m_int_hat = m_int_hat, sim_est = sim_est))
#     }) 
#     
#   })
#   names(tree_est_list) <- c('Whittle', 'periodogram') 
#   #saveRDS(object = tree_est_list, file = "C:\\Users\\howard\\Desktop\\PhD Research\\LRD Resampling\\Proj 2 Estimation of Hermite rank\\Code\\cache\\tree_est_list.rds")
#   
#   
#   tree_est_list <- readRDS(file = "cache\\tree_est_list.rds")
#   ### Hermite Hermite Rank Estimation Comparison Plot for Xt vs Xt^2
#   plot_names <- c('Samplewise Comparison (Local Whittle)',
#                   'Samplewise Comparison (Log-periodogram)')
#   hist_names <- c('Histogram Comparison for Tree Rings (Local Whittle)',
#                   'Histogram Comparison for Tree Rings (Log-periodogram)')
#   sim_names <- c('Histogram Comparison for Simulations (Local Whittle)',
#                  'Histogram Comparison for Simulations (Log-periodogram)')
#   
#   ### for two different am estimation methods: Whittle, log-periodogram
#   Treering_plot_list <- list()
#   tree_sim_comparison <- list()
#   for(i in c(1,2)){
#     tree_sim_comparison[[i]] <- list()
#     ### count the frequency of m_hat(Xt) > m_hat(Xt^2)
#     tree_sim_comparison[[i]]$tree_m_comparison <- sapply(tree_est_list[[i]], function(lst){
#       lst$m_int_hat
#     }) %>% 
#       t() %>% 
#       as_tibble() %>% 
#       `colnames<-`(c('Xt', 'Xt2')) %>% 
#       apply(1, function(x){
#         if(x[2]>x[1]){
#           return(1)
#         }else if(x[2] == x[1]){
#           return(0)
#         }else{
#           return(-1)
#         }
#       }) %>% 
#       table()
#     
#     # m_comparision for simulations
#     tree_sim_comparison[[i]]$sim_m_comparison <- sapply(tree_est_list[[i]], function(lst){
#       lst$sim_est %>% 
#         apply(1, function(x){
#           if(x[2]>x[1]){
#             return(1)
#           }else if(x[2] == x[1]){
#             return(0)
#           }else{
#             return(-1)
#           }
#         }) %>% 
#         table()
#     }) %>% 
#       t() %>% 
#       colSums()
#     
#     ### plots for Hermite rank estimation of tree rings data
#     # remark: type 1 for Xt and type 2 for Xt^2
#     tree_m_hat_tbl <- sapply(tree_est_list[[i]], function(lst){
#       lst$m_int_hat
#     }) %>% 
#       t() %>% 
#       as_tibble() %>% 
#       `colnames<-`(c(1,2)) %>% 
#       pivot_longer(cols = everything(), names_to = 'type', values_to = 'm_hat') %>% 
#       dplyr::mutate(m_hat = ifelse(m_hat>20, 20, m_hat), 
#                     index = rep.int(rep(seq(1, 0.5*n())), 
#                                         times = rep(2, times= 0.5*n())))
#     
#     Treering_plot_list[[i]] <- list()
#     Treering_plot_list[[i]][[1]] <- ggplot(tree_m_hat_tbl, aes(x = index, y = m_hat, colour = type)) +
#       geom_point() + 
#       scale_colour_manual(name = '',
#                          values=c("#69b3a2", "#404080"), 
#                          labels= lapply(c('$X_t$', '$X_t^2$'), TeX)) + 
#       scale_y_continuous(breaks = seq(0, 20, by = 5), lim = c(0, 20), label = c(seq(0, 15, by = 5), '20+')) +
#       ylab(TeX('$\\widehat{m}$')) +
#       ggtitle(plot_names[i]) + 
#       theme(plot.title = element_text(hjust = 0.5)) 
#     
#     tree_m_hat_tbl <- tree_m_hat_tbl %>% 
#       dplyr::mutate(m_hat_trunc = ifelse(m_hat>=10, 10, m_hat))
#     
#     Treering_plot_list[[i]][[2]] <- ggplot(tree_m_hat_tbl , aes(x=m_hat_trunc, fill=type)) +
#       geom_histogram( color="#e9ecef", alpha=0.4, position = 'identity', binwidth = 1) +
#       guides(fill=guide_legend(title=NULL)) + 
#       scale_fill_manual(values=c("#69b3a2", "#404080"), labels=lapply(c('$X_t$', '$X_t^2$'), TeX)) +
#       scale_x_continuous(breaks = seq(1, 10), lim = c(0, 11), labels = c(1:9, '9+')) +
#       scale_y_continuous(breaks = c(0, 100, 200, 300), lim = c(0, 300)) +
#       xlab(TeX('$\\widehat{m}$')) +
#       ggtitle(hist_names[i]) + 
#       theme(plot.title = element_text(hjust = 0.5))
#     
#     ### plots for Hermite rank estimation of simulated FGuassian
#     sim_m_hat_tbl <- lapply(tree_est_list[[i]], function(lst){
#       lst$sim_est
#     }) %>% 
#       Reduce(f = rbind) %>% 
#       apply(2, function(m){
#         ifelse(m>=10, 10, m)
#       }) %>% 
#       as_tibble() %>% 
#       `colnames<-`(c(1,2)) %>% 
#       pivot_longer(cols = everything(), names_to = 'type', values_to = 'm_hat')
#     
#     Treering_plot_list[[i]][[3]] <- ggplot(sim_m_hat_tbl , aes(x=m_hat, fill=type)) +
#       geom_histogram( color="#e9ecef", alpha=0.4, position = 'identity', binwidth = 1) +
#       guides(fill=guide_legend(title=NULL)) + 
#       scale_fill_manual(values=c("#69b3a2", "#404080"), labels=lapply(c('$X_t$', '$X_t^2$'), TeX)) +
#       scale_x_continuous(breaks = seq(1, 10), lim = c(0, 11), labels = c(1:9, '9+')) +
#       scale_y_continuous(breaks = c(0, 10000, 20000, 30000), lim = c(0, 32000)) +
#       xlab(TeX('$\\widehat{m}$')) +
#       ggtitle(sim_names[i]) + 
#       theme(plot.title = element_text(hjust = 0.5))
#     
#     ### Chi-square test for comparing two histograms
#     sim_m_hat_freq <- sim_m_hat_tbl %>% 
#       group_by(type, m_hat) %>% 
#       summarise(n = n()) %>% 
#       mutate(freq = n/sum(n))
#     
#     tree_m_hat_freq <- tree_m_hat_tbl %>% 
#       group_by(type, m_hat_trunc) %>% 
#       summarise(n = n()) %>% 
#       mutate(freq = n/sum(n))
#     
#     tree_sim_comparison[[i]]$Chisq_pvalues <- sapply(c(1, 2), function(j){
#       sim_subset <- sim_m_hat_freq %>% 
#         filter(type == j) %>% 
#         select(freq, n)
#       tree_subset <- tree_m_hat_freq %>% 
#         filter(type == j) %>% 
#         select(freq, n)
#       
#       # chi-square test between hist of Simulations vs tree ring for Xt and Xt^2
#       chi_stat <- sum((sim_subset$freq - tree_subset$freq)^2/
#                         (sim_subset$freq/sum(sim_subset$n) + tree_subset$freq/sum(tree_subset$n)))
#       p_value <- 1- pchisq(chi_stat, df = 9, lower.tail = T)
#     }) %>% 
#       `names<-`(c('Xt','Xt2'))
#   }
#   
#   #tree_sim_comparison
#   ### merge plots
#   TreeRing_plot <- grid.arrange(Treering_plot_list[[1]][[1]], Treering_plot_list[[2]][[1]],
#                                 Treering_plot_list[[1]][[2]], Treering_plot_list[[2]][[2]],
#                                 Treering_plot_list[[1]][[3]], Treering_plot_list[[2]][[3]], nrow = 3)
#   
#   ggsave("TreeRing_plot.png", plot = TreeRing_plot,  width = 11, height = 10.5)
# }




