library(EQL)
library(nortest)
library(zoo)
library(LongMemoryTS)
library(gridExtra)
library(grid)
library(latex2exp)

source("functions/parallel.R")
source("functions/Hermite_rank_estimator.R")
source("functions/Subsampling_mbb_Hermite_test.R")

### Hermite rank test based on Hermite rank estimator
hermite_test_byest <-
  function(x, l, k_array, am_option, m0, B, alpha_sig) {
    ####
    # x: a single time series
    # l: block size
    # k_array: different k for squared blocks
    # am_option: 'Whittle' or 'periodogram'
    # m0: for H0: m=m0
    # B: number of bootstrap samples for CI
    # alpha_sig: significant level for the test
    ####

    # estimate alpha*m
    n <- length(x)
    am_est_func <- function(x, n) {
      if (am_option == 'Whittle') {
        am_hat <- (1 - 2 * GSE(x,  m = floor(1 + n ^ 0.65)))
      } else if (am_option == 'periodogram') {
        am_hat <- (1 - 2 * gph(x,  m = floor(1 + n ^ 0.7)))
      }
      return(am_hat)
    }
    am_hat <- am_est_func(x = x, n = n)
    am_hat <- max(0.01, min(0.99, abs(am_hat)))
    a_hat <- max(0.01, min(0.99, am_hat / m0))
    m_hat <-
      HermiteRank_est_func(
        x = x,
        l = l,
        k_array = k_array,
        am_option = am_option
      )$m_hat
    
    # simulate B samples given a_hat
    m_hat_stars_CI <-
      LRD_gen(
        n = n,
        alpha = a_hat,
        nrep = B,
        option = 'FGaussian'
      ) %>%
      sapply(function(x) {
        # browser()
        xt_star <- hermite(x = x, n = m0, prob = T)
        #
        am_star <- am_est_func(x = xt_star, n = n)
        am_star <- max(0.01, min(0.99, abs(am_star)))
        m_hat_star <- am_star/a_hat
        
        return(m_hat_star)
      }) %>%
      quantile(probs = c(0, 1 - alpha_sig))
    
    # decision = T, reject H0
    decision <- (m_hat >= m_hat_stars_CI[2])
    
    return(decision)
  }


hermite_test_simulation <-
  function(n,
           alpha,
           G_list,
           nrep,
           m0,
           B,
           cores,
           alpha_sig = 0.05) {
    ####
    # n: length for the simulated ts
    # alpha: long-memory exponent
    # G_list: a list of transformation G(.)
    # nrep: number of replications
    # k_array: different k for squared blocks
    # am_option: 'Whittle' or 'periodogram'
    # m0: for H0: m=m0
    ####
    
    zt_list <- LRD_gen(
      n = n,
      alpha = alpha,
      nrep = nrep,
      option = 'FGaussian'
    )
   
    # parameters for Hermite estimation
    l <- round(2 * n ^ (1 / 4))
    k_max <- 3 * (n / l) ^ (1 / 5)
    k_min <- max(2, round(k_max ^ (1 / 2)))
    k_array <- seq(k_min, max(k_min + 1, round(k_max)))
    
    
    nrep_bools <- mclapply(zt_list, function(zt) {
      test_bools <- sapply(seq(1, length(G_list)), function(i) {
        G <- G_list[[i]]
        ts <- G(zt)
        test_bool <-
          hermite_test_byest(
            x = ts,
            l = l,
            k_array = k_array,
            am_option = 'Whittle',
            m0 = m0,
            B = B,
            alpha_sig = alpha_sig
          )

        return(test_bool)
      })

      return(test_bools)
    }, cores=cores) %>%
      Reduce(f = rbind)
    
    reject_prob <- colMeans(nrep_bools)
    
    return(reject_prob)
  }


if (sys.nframe() == 0) {
  G_list <- list()
  G_list[[1]] <- function(x) sin(x)
  G_list[[2]] <- function(x)
    hermite(x = x, n = 2, prob = T)  + hermite(x = x, n = 4, prob = T) + + hermite(x = x, n = 5, prob = T)
  G_list[[3]] <- function(x)
      hermite(x = x, n = 3, prob = T)  + hermite(x = x, n = 3, prob = T)
  
  n0_array <- c(800, 2000, 5000)
  alpha_array <- c(0.1)
  B_array <- 50
  nrep <- 500
  
  # JASA method
  set.seed(11)
  summary_all_hr <- lapply(B_array, function(B) {
    summary_all_B <- lapply(alpha_array, function(alpha) {
      summary_n0s <- lapply(n0_array, function(n) {
        print(paste0('B=', B, ', a=', alpha, ', n=', n))
        rej_probs <-
          hermite_test_simulation(
            n = n,
            alpha = alpha,
            G_list = G_list,
            nrep = nrep,
            m0 = 1,
            B = B,
            cores = max(2, detectCores() - 6),
            alpha_sig = 0.05
          )
        summary_vec <- c(alpha, n, rej_probs)
        names(summary_vec) <- c('alpha', 'n', 'm=1', 'm=2', 'm=3')
        
        return(summary_vec)
      }) %>%
        Reduce(f = rbind)
    }) %>%
      Reduce(f = rbind)
    
    return(summary_all_B)
  })
  
  # JCGS method
  summary_all_jcgs <- lapply(B_array, function(B){
    summary_all_B <- lapply(alpha_array, function(alpha){
      summary_n0s <- lapply(n0_array, function(n){
        rej_probs <- two_sample_test_simulation(n=n, alpha=alpha, l=floor(n^(1/2)), B=B, nrep=nrep, G_list=G_list, 
                                                cores = max(2, detectCores()-6))
        summary_vec <- c(alpha, n, rej_probs)
        names(summary_vec) <- c('alpha', 'n', 'm=1', 'm=2', 'm=3')
        
        return(summary_vec)
      }) %>% 
        Reduce(f=rbind)
    }) %>% 
      Reduce(f=rbind)
    
    return(summary_all_B)
  })
  
  summary_all_merged <- rbind(summary_all_hr[[1]], summary_all_jcgs[[1]]) %>% 
    as_tibble() %>% 
    mutate('type' = rep(c('hr', 'jcgs'), each =3)) 
  
  # saveRDS(summary_all_merged, file="cache/test_comparison_rejections_seed11.rds")
  # saveRDS(summary_all_merged, file="cache/test_comparison_rejections.rds")
}



####################### Hermite rank estimator diagnostic
hermite_rank_diag <- function(ts, m0_array, by, B, cores){
  ####
  # ts: a vector of time series observations
  # m0_array: a vector of m0 for diagnostic
  # by: int, increment of k_array
  # B: int, number of bootstrap sample
  # cores: number of cores for parallel computing
  ####
  
  n <- length(ts)
  
  # parameters for Hermite estimation
  l <- round(2 * n ^ (1 / 4))
  k_max <- round(3 * (n / l) ^ (1 / 5))
  k_min <- max(2, round(k_max ^ (1 / 2))+1)
  k_array <- seq(k_min, max(k_min + 1, round(k_max)))
  cat('k_min =',k_min, ', k_max =',k_max, '\n')
  est_df <- HermiteRank_est_func(x = ts, l = l, 
                                k_array = seq(k_min, k_max, by = by),
                                am_option = 'Whittle')

  m_hat <- est_df$m_hat
  am_hat <- min(0.99, max(0.01, abs(est_df$am_hat)))
  a_hat <- max(0.01, min(0.99, am_hat / m_hat))
  
  xt_summary <- tibble('m_hat'=m_hat, 'a_hat'=a_hat)
  
  m_stars_all <- mclapply(m0_array, function(m0){
    a_hat0 <- am_hat/m0
    
    m_stars_m0 <- LRD_gen(
      n = n,
      alpha = a_hat0,
      nrep = B,
      option = 'FGaussian'
    ) %>% 
      sapply(function(z){
        xt_star <- hermite(x = z, n = m0, prob = T) 
        am_star <- (1 - 2 * GSE(X=xt_star,  m = floor(1 + n ^ 0.65)))
        am_star <- min(0.99, max(0.01, abs(am_star)))
        m_star <- HermiteRank_est_func(x=xt_star, l=l, k_array=k_array, am_option = am_option)$m_hat
        a_star <- am_star/m_star
        
        m_stars <- c()
        m_stars[1] <- m0
        m_stars[2] <- max(1, min(am_star/a_hat0, 15))
        
        return(m_stars)
      }) %>% 
      t()
    
    return(m_stars_m0)
  }, cores=cores) %>% 
    Reduce(f=rbind) %>% 
    as_tibble() %>% 
    `colnames<-`(c('m0','m_star1'))
  
  return(list(xt_summary = xt_summary, m_stars_all = m_stars_all))
}

####### DJI diagnostic

if(sys.nframe()==0){

  m0_array <- c(1,2,3, 5, 7, 10)
  B <- 200
  
  set.seed(11)
  m_diag_DJI <- hermite_rank_diag(ts=DJI$X_t, m0_array=m0_array, by=2, B=B, cores=10)  
  
  m_diag_DJI$m_stars_all$m0_parsed <- factor(
    m_diag_DJI$m_stars_all$m0,
    labels = c(
      '1' = parse(text = TeX('m_0=1')),
      '2' = parse(text = TeX('m_0=2')),
      '3' = parse(text = TeX('m_0=3')),
      '5' = parse(text = TeX('m_0=5')),
      '7' = parse(text = TeX('m_0=7')),
      '10' = parse(text = TeX('m_0=10'))
    )
  )
  
  m_diag_DJI_hists <- m_diag_DJI$m_stars_all %>%
    ggplot(aes(x = m_star1)) +
    geom_histogram(aes(y = ..density..), color='black', fill=c("#69b3a2"), binwidth = 1) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      color = "red",
      size = 1
    ) +
    ylim(0, 1) +
    ylab('Frequency') + 
    xlab(TeX('$\\widehat{m}^*$')) +
    scale_x_continuous(breaks = seq(1, 15, by = 2), lim = c(0, 15), label = c(seq(1, 14, by = 2), '15+')) + 
    ggtitle('Histograms of Bootstrap Estimations (Dow Jones Index)') +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap( ~ m0_parsed, labeller=label_parsed, scales = "free"
    ) + 
    theme(axis.title = element_text(size = 12),
          strip.text = element_text(size = 12))
  
  # m_diag_SnP500 <- hermite_rank_diag(ts=SnP500$X_t, m0_array=m0_array, by=1, B=B, cores=10)    
  # m_diag_SnP500$m_stars_all$m0_parsed <- factor(
  #   m_diag_SnP500$m_stars_all$m0,
  #   labels = c(
  #     '1' = parse(text = TeX('m_0=1')),
  #     '2' = parse(text = TeX('m_0=2')),
  #     '3' = parse(text = TeX('m_0=3')),
  #     '5' = parse(text = TeX('m_0=5')),
  #     '7' = parse(text = TeX('m_0=7')),
  #     '10' = parse(text = TeX('m_0=10'))
  #   )
  # )
  # 
  # m_diag_SnP500_hists <- m_diag_SnP500$m_stars_all %>%
  #   ggplot(aes(x = m_star1)) +
  #   geom_histogram(aes(y = ..density..), color='black', fill=c("#69b3a2"), binwidth = 1) +
  #   geom_vline(
  #     xintercept = 1,
  #     linetype = "dashed",
  #     color = "red",
  #     size = 1
  #   ) +
  #   ylim(0, 1) +
  #   ylab('Frequency') + 
  #   xlab(TeX('$\\widehat{m}^*$')) +
  #   scale_x_continuous(breaks = seq(1, 15, by = 2), lim = c(0, 15), label = c(seq(1, 14, by = 2), '15+')) + 
  #   ggtitle('Histograms of Bootstrap Estimations (S&P 500)') +
  #   theme(plot.title = element_text(hjust = 0.5)) +
  #   facet_wrap( ~ m0_parsed, labeller=label_parsed, scales = "free"
  #   ) + 
  #   theme(axis.title = element_text(size = 12),
  #         strip.text = element_text(size = 12))

  ggsave("plots/DJI_diagnostic_histograms.png", plot = m_diag_DJI_hists,  width = 11, height = 7)
  # ggsave("plots/SnP500_diagnostic_histograms.png", plot = m_diag_SnP500_hists,  width = 11, height = 7)
}


























# 
# #####################################################
# if(sys.nframe()==0){
#   G_list <- list()
#   G_list[[1]] <- function(x) sin(x)
#   G_list[[2]] <-
#     function(x)
#       hermite(x = x, n = 2, prob = T)  + hermite(x = x, n = 3, prob = T)
#   
#   n_array <- c(1000, 5000)
#   B <- 100
#   m0_array <- seq(1, 9, by=1)
#   
#   m_diag_list_seeds <- lapply(c(1:50), function(i){
#     set.seed(i)
#     m_diag_list <- lapply(G_list, function(G){
#       m_diag_n <- lapply(n_array, function(n){
#         m_diag <- hermite_rank_diag(n=n, alpha=0.4, G=G, m0_array=m0_array, B=B, am_option='Whittle', cores=10)     
#       })
#     })
#   })
# 
#   for(i in c(1:50)){
#     p1 <- m_diag_list_seeds[[i]][[1]][[1]] %>%
#       Reduce(f=rbind) %>%
#       t() %>%
#       as_tibble() %>%
#       `colnames<-`(paste0('m0=', m0_array)) %>%
#       pivot_longer(everything(),
#                    names_to = "m0",
#                    values_to = "m_est") %>%
#       ggplot(aes(x=m_est, fill = m0)) +
#       geom_histogram(binwidth = 1, aes(y = ..density..)) +
#       ylim(0, 1) + 
#       ggtitle('G1, n=1000') + 
#       theme(plot.title = element_text(hjust = 0.5)) +
#       facet_wrap(~m0)
#     
#     p2 <- m_diag_list_seeds[[i]][[1]][[2]] %>%
#       Reduce(f=rbind) %>%
#       t() %>%
#       as_tibble() %>%
#       `colnames<-`(paste0('m0=', m0_array)) %>%
#       pivot_longer(everything(),
#                    names_to = "m0",
#                    values_to = "m_est") %>%
#       ggplot(aes(x=m_est, fill = m0)) +
#       geom_histogram(binwidth = 1, aes(y = ..density..)) +
#       ylim(0, 1) + 
#       ggtitle('G1, n=5000') + 
#       theme(plot.title = element_text(hjust = 0.5)) +
#       facet_wrap(~m0)
#     
#     p3 <- m_diag_list_seeds[[i]][[2]][[1]] %>%
#       Reduce(f=rbind) %>%
#       t() %>%
#       as_tibble() %>%
#       `colnames<-`(paste0('m0=', m0_array)) %>%
#       pivot_longer(everything(),
#                    names_to = "m0",
#                    values_to = "m_est") %>%
#       ggplot(aes(x=m_est, fill = m0)) +
#       geom_histogram(binwidth = 1, aes(y = ..density..)) +
#       ylim(0, 1) + 
#       ggtitle('G2, n=1000') + 
#       theme(plot.title = element_text(hjust = 0.5)) +
#       facet_wrap(~m0)
#     
#     p4 <- m_diag_list_seeds[[i]][[2]][[2]] %>%
#       Reduce(f=rbind) %>%
#       t() %>%
#       as_tibble() %>%
#       `colnames<-`(paste0('m0=', m0_array)) %>%
#       pivot_longer(everything(),
#                    names_to = "m0",
#                    values_to = "m_est") %>%
#       ggplot(aes(x=m_est, fill = m0)) +
#       geom_histogram(binwidth = 1, aes(y = ..density..)) +
#       ylim(0, 1) + 
#       ggtitle('G2, n=5000') + 
#       theme(plot.title = element_text(hjust = 0.5)) +
#       facet_wrap(~m0)
#     
#     p_combined <- grid.arrange(p1, p2, p3, p4, ncol = 2)
#     
#     ggsave(paste0("plots/temp/diagnostic_", i,".png"), plot = p_combined,  width = 11*2, height = 7*2)
#   }
# 
#   # m_diag %>% 
#   #   Reduce(f=rbind) %>% 
#   #   t() %>% 
#   #   as_tibble() %>% 
#   #   `colnames<-`(c(1:10)) %>% 
#   #   pivot_longer(everything(),
#   #               names_to = "m0",
#   #               values_to = "m_est") %>% 
#   #   ggplot(aes(x=m_est)) +
#   #   geom_histogram(alpha=0.6, binwidth = 1) +
#   #   facet_wrap(~m0)
#   # 
# 
#  
# }
# 
# 
# 
# 
# 
# if(sys.nframe()==0){
#   G_list <- list()
#   G_list[[1]] <- function(x) hermite(x = x, n = 1, prob = T)  + hermite(x = x, n = 2, prob = T)
#   G_list[[2]] <-
#     function(x)
#       hermite(x = x, n = 2, prob = T) + hermite(x = x, n = 3, prob = T) + hermite(x = x, n = 4, prob = T)
#   
#   n_array <- c(1000, 5000)
#   B <- 100
#   m0_array <- seq(1, 9, by=1)
#   
#   m_diag_list_seeds_2 <- lapply(c(1:50), function(i){
#     set.seed(i)
#     m_diag_list <- lapply(G_list, function(G){
#       m_diag_n <- lapply(n_array, function(n){
#         m_diag <- hermite_rank_diag(n=n, alpha=0.45, G=G, m0_array=m0_array, B=B, am_option='Whittle', cores=10)     
#       })
#     })
#   })
#   
#   for(i in c(1:50)){
#     p1 <- m_diag_list_seeds_2[[i]][[1]][[1]] %>%
#       Reduce(f=rbind) %>%
#       t() %>%
#       as_tibble() %>%
#       `colnames<-`(paste0('m0=', m0_array)) %>%
#       pivot_longer(everything(),
#                    names_to = "m0",
#                    values_to = "m_est") %>%
#       ggplot(aes(x=m_est, fill = m0)) +
#       geom_histogram(binwidth = 1, aes(y = ..density..)) +
#       ylim(0, 1) + 
#       ggtitle('G1, n=1000') + 
#       theme(plot.title = element_text(hjust = 0.5)) +
#       facet_wrap(~m0)
#     
#     p2 <- m_diag_list_seeds_2[[i]][[1]][[2]] %>%
#       Reduce(f=rbind) %>%
#       t() %>%
#       as_tibble() %>%
#       `colnames<-`(paste0('m0=', m0_array)) %>%
#       pivot_longer(everything(),
#                    names_to = "m0",
#                    values_to = "m_est") %>%
#       ggplot(aes(x=m_est, fill = m0)) +
#       geom_histogram(binwidth = 1, aes(y = ..density..)) +
#       ylim(0, 1) + 
#       ggtitle('G1, n=5000') + 
#       theme(plot.title = element_text(hjust = 0.5)) +
#       facet_wrap(~m0)
#     
#     p3 <- m_diag_list_seeds_2[[i]][[2]][[1]] %>%
#       Reduce(f=rbind) %>%
#       t() %>%
#       as_tibble() %>%
#       `colnames<-`(paste0('m0=', m0_array)) %>%
#       pivot_longer(everything(),
#                    names_to = "m0",
#                    values_to = "m_est") %>%
#       ggplot(aes(x=m_est, fill = m0)) +
#       geom_histogram(binwidth = 1, aes(y = ..density..)) +
#       ylim(0, 1) + 
#       ggtitle('G2, n=1000') + 
#       theme(plot.title = element_text(hjust = 0.5)) +
#       facet_wrap(~m0)
#     
#     p4 <- m_diag_list_seeds_2[[i]][[2]][[2]] %>%
#       Reduce(f=rbind) %>%
#       t() %>%
#       as_tibble() %>%
#       `colnames<-`(paste0('m0=', m0_array)) %>%
#       pivot_longer(everything(),
#                    names_to = "m0",
#                    values_to = "m_est") %>%
#       ggplot(aes(x=m_est, fill = m0)) +
#       geom_histogram(binwidth = 1, aes(y = ..density..)) +
#       ylim(0, 1) + 
#       ggtitle('G2, n=5000') + 
#       theme(plot.title = element_text(hjust = 0.5)) +
#       facet_wrap(~m0)
#     
#     p_combined <- grid.arrange(p1, p2, p3, p4, ncol = 2)
#     
#     ggsave(paste0("plots/temp/diagnostic_new_", i,".png"), plot = p_combined,  width = 11*2, height = 7*2)
#   }
#   
#   
# }
# 
# 
