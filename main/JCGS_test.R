library(tidyverse)
library(EQL)
library(nortest)
library(zoo)
library(LongMemoryTS)
library(gridExtra)
library(grid)

source("functions\\parallel.R")
source("functions\\LRD_Generator_new.R")

# generate B standardized mean based on moving block bootstrap samples 
mvblock_bt <- function(ts, l, am,  B){
  ###########
  # ts: a vector
  # l: block size for moving block bootstrap
  # am: long-memory parameter alpha*m
  # B: number of bootstrap samples
  ###########
  n <- length(ts)
  k <- floor(n/l)
  block_means <- rollmean(x=ts, k=l)
  mu_hat <- mean(block_means)  
  
  # sample indices for B bootstrap samples
  indices <- sample(x=seq(1, n-l+1), size=B*k, replace=T)
  bt_stats <- sapply(seq(0, B-1), function(i){
    ind <- indices[seq(k*i+1, k*i+k)]
    bt_stat <- k^(1/2)*l^(am/2)*(mean(block_means[ind]) - mu_hat)
    
    return(bt_stat)
  })
  return(bt_stats)
}

# generate B standardized mean via subsampling method
subsample <- function(ts, l, am, B){
  ###########
  # ts: a vector
  # l: block size for moving block bootstrap
  # am: long-memory parameter alpha*m
  # B: number of bootstrap samples
  ###########
  n <- length(ts)
  k <- floor(n/l)
  block_means <- rollmean(x=ts, k=l)
  mu_hat <- mean(block_means)  
  
  # select B blocks
  indices <- sample(x=seq(1, n-l+1), size=B, replace=T)
  sub_stats <- sapply(seq(1, B), function(i){
    sub_stat <- l^(am/2)*(block_means[indices[i]] - mu_hat)
    
    return(sub_stat)
  })
  return(sub_stats)
}

# Hermite rank test for m=1 by two sample test (subsample vs MBB)
hermite_two_sample_test <- function(ts, l, B, alpha_sig, am_option = 'Whittle'){
  ###########
  # ts: a vector
  # l: block size for moving block bootstrap
  # B: number of bootstrap samples
  # alpha_sig: significant level alpha=0.05
  # am_option: 'Whittle' or 'periodogram'
  ########### 
  
  assert_that(am_option %in% c('Whittle', 'periodogram'))
  
  # estimates for alpha*m
  am_est_func <- function(x, n){
    if(am_option == 'Whittle'){
      am_hat <- (1 - 2 *GSE(x,  m=floor(1+n^0.65)))
    }else if(am_option == 'periodogram'){
      am_hat <- (1 - 2 *gph(x,  m=floor(1+n^0.7)))
    }
    am_hat <- max(0.01, min(0.99, abs(am_hat)))
    return(am_hat)
  }
  
  am_hat <- am_est_func(x=ts, n=length(ts))
  bt_stats <- mvblock_bt(ts=ts, l=l, am=am_hat, B=B)
  sub_stats <- subsample(ts=ts, l=l, am=am_hat, B=B)
  
  ks_pvalue <- ks.test(bt_stats, sub_stats)$p.value
  return(ks_pvalue<alpha_sig)
}

# Hermite rank test simulator for computing rejected probs
two_sample_test_simulation <- function(n, alpha, l, B, nrep, G_list, option = 'FGaussian', am_option = 'Whittle', alpha_sig=0.05, cores=2){
  ###########
  # n: lenth of generated ts
  # ts: a vector
  # alpha: memory exponent for Fractional Gaussian Zt
  # l: block size for moving block bootstrap
  # B: number of bootstrap samples
  # nrep: number of generates Xt = G(Zt)
  # G_list: a list of functions G(.)
  # alpha_sig: significant level alpha=0.05
  # am_option: 'Whittle' or 'periodogram'
  # option: 'FGaussian' or 'Farima'
  # cores: number of cores for parallel computing
  ########### 
  
  ts_list <- LRD_gen(
    n = n,
    alpha = alpha,
    nrep = nrep,
    option = option
  )
  
  nrep_bools <- mclapply(ts_list, function(zt){
    test_bools <- sapply(G_list, function(G){
      ts <- G(zt)
      test_bool <- hermite_two_sample_test(ts=ts, l=l, B=B, alpha_sig=alpha_sig, am_option=am_option)
      
      return(test_bool)
    })
    return(test_bools)
  }, cores=cores) %>% 
    Reduce(f = rbind)
  
  reject_prob <- colMeans(nrep_bools)
  return(reject_prob)
}

# function for coercing multiple legends to one legend for grid.arrange
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(0.98, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  # invisible(combined)
  return(combined)
}


if(sys.nframe()==0){
  G_list <- list()
  G_list[[1]] <- sin
  G_list[[2]] <- function(x) hermite(x=x, n=2, prob = T) + hermite(x=x, n=4, prob = T) + hermite(x=x, n=5, prob = T)
  G_list[[3]] <- function(x) hermite(x=x, n=3, prob = T) + hermite(x=x, n=4, prob = T)
  
  n0_array <- c(800, 2000, 5000)
  alpha_array <- c(0.15, 0.3)
  B_array <- seq(50, 500, by=50)
  nrep <- 500
  
  # compute test rejected probs over combinations of B, alpha, n
  set.seed(44)
  summary_all <- lapply(B_array, function(B){
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

  rej_probs_df <- summary_all %>% 
    Reduce(f=rbind) %>% 
    cbind('B'=rep(B_array, each = 6)) %>% 
    as_tibble() %>% 
    pivot_longer(!c('alpha', 'n', 'B'), names_to = "m", values_to = "prob")
  
  # saveRDS(rej_probs_df, 'cache/rej_probs_df.rds')
  
  # draw plots for the paper
  plot_list <- list()
  for(i in c(1:3)){
    n0 <- n0_array[i]
    plot_list[[i]] <- rej_probs_df %>%
      filter(n==n0) %>%
      ggplot(aes(x=B, y=prob, linetype=factor(alpha), shape=m,
                 group=interaction(alpha, m))) +
      geom_line()+
      geom_point(size = 2) +
      guides(linetype=guide_legend(title='Long-memory', label.position = 'right', order = 1),
             shape=guide_legend(title='Rank', order = 2)) +
      xlim(c(0, 500)) +
      ylim(c(0, 0.7)) +
      xlab(TeX('B')) +
      ylab(TeX('Rejection Probablity')) +
      ggtitle(paste('Hermite Rank Test (n=', n0, ')', sep = '')) +
      theme(plot.title = element_text(hjust = 0.5))  +
      scale_linetype_discrete(labels=lapply(sprintf('$\\alpha = %s$',
                                                    round(alpha_array, digits = 3)), TeX)) + 
      scale_shape_discrete(labels=lapply(sprintf('$m = %s$',
                                                  round(c(1,2,3), digits = 1)), TeX)) +
      geom_hline(yintercept=0.05, linetype="dashed", 
                 color = "red", size=1)
  }
  
  B0 <- 50
  plot_list[[4]] <- rej_probs_df %>%
    filter(B==B0) %>%
    ggplot(aes(x=factor(n), y=prob, linetype=factor(alpha), shape=m,
               group=interaction(alpha, m))) +
    geom_line()+
    geom_point(size = 2) +
    guides(linetype=guide_legend(title='Memory', label.position = 'right', order = 1),
           shape=guide_legend(title='Rank', order = 2)) +
    # ylim(c(0, 0.25)) +
    xlab(TeX('n')) +
    ylab(TeX('Rejection Probablity')) +
    ggtitle(paste('Hermite Rank Test (B=', B0, ')', sep = '')) +
    theme(plot.title = element_text(hjust = 0.5))  +
    scale_linetype_discrete(labels=lapply(sprintf('$\\alpha=%s$',
                                                  round(alpha_array, digits = 3)), TeX)) + 
    geom_hline(yintercept=0.05, linetype="dashed", 
               color = "red", size=1)

  p <- grid_arrange_shared_legend(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], ncol = 2, nrow = 2, position='right')
  ggsave("plots/Hermite_two_sample_test.png", plot = p,  width = 11, height = 7)
  
}













