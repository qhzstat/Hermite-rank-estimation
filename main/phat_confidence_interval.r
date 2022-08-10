library(tidyverse)
library(latex2exp)
library(gridExtra)
library(grid)

# load tree_est data
tree_est_list <- readRDS(file = "cache\\tree_est_list.rds")
tree_ranks_list <- tree_est_list$Whittle

if(sys.nframe()==0){
  # Histogram comparison between \hat{m}^{(2)} and mixture of Z_t and Z_t^2 
  # select the results by Whittle estimator
  tree_m_hat_tbl <- sapply(tree_ranks_list, function(lst){
    lst$m_int_hat
  }) %>% 
    t() %>% 
    as_tibble() %>% 
    `colnames<-`(c(1,2)) %>% 
    pivot_longer(cols = everything(), names_to = 'type', values_to = 'm_hat') %>% 
    dplyr::mutate(m_hat = ifelse(m_hat>20, 20, m_hat), 
                  index = rep.int(rep(seq(1, 0.5*n())), 
                                  times = rep(2, times= 0.5*n()))) %>% 
    mutate(m_hat_trunc = ifelse(m_hat>=10, 10, m_hat))
  
  
  sim_m_hat_tbl <- lapply(tree_ranks_list, function(lst){
    lst$sim_est
  }) %>% 
    Reduce(f = rbind) %>% 
    apply(2, function(m){
      ifelse(m>=10, 10, m)
    }) %>% 
    as_tibble() %>% 
    `colnames<-`(c(1,2)) %>% 
    pivot_longer(cols = everything(), names_to = 'type', values_to = 'm_hat')
  
  
  # Tree ring's table for the proportion of m_hat for X_t^2
  tree_xt2 <- tree_m_hat_tbl %>%
    filter(type == '2') %>%
    group_by(m_hat_trunc) %>%
    summarise(count = n()) %>%
    ungroup() %>%
    summarise(m_hat = m_hat_trunc,
              count = count,
              proportion = count / sum(count),
              source = 'Tree Ring')
  
  # Contrast group's table for the proportion of m_hat for the mixtures (w*X_t + (1-w)*X_t^2)
  weight <- 0.5
  sim_wavg <- sim_m_hat_tbl %>%
    group_by(type, m_hat) %>%
    summarise(count = n()) %>%
    ungroup(m_hat) %>%
    summarise(m_hat = m_hat,
              proportion = count / sum(count),
              count = count) %>%
    group_by(m_hat) %>%
    summarise(
      count = count[1],
      proportion = weight * proportion[1] + (1 - weight) * proportion[2],
      source = 'Simulation'
    ) 
  
  squared_tree_ring_plot <- rbind(tree_xt2, sim_wavg) %>%
    ggplot(aes(x = m_hat, y = proportion, fill = source)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_manual(
      name = 'Source',
      values = c("#69b3a2", "#404080"),
      labels = lapply(c('$Mixture$', '$Tree\\, Ring$'), TeX)
    ) +
    scale_x_continuous(
      breaks = seq(1, 10),
      lim = c(0, 11),
      labels = c(1:9, '9+')
    ) +
    xlab(TeX('$\\widehat{m}$')) +
    ylab('Density') +
    ylim(c(0, 0.6)) +
    ggtitle(
      TeX(
        '$Hermite\\, Rank\\, of\\, Squared\\, Tree\\, Rings\\, X_t^2\\, vs\\, Mixture$'
      )
    ) +
    theme(plot.title = element_text(hjust = 0.5))

  # chi-square test between hist of m_hat for the mixtures (w*X_t + (1-w)*X_t^2)
  chi_stat <- sum((sim_wavg$proportion - tree_xt2$proportion)^2/
                    (sim_wavg$proportion/sum(sim_wavg$count) + tree_xt2$proportion/sum(tree_xt2$count)))
  p_value <- 1- pchisq(chi_stat, df = 9, lower.tail = T)
  
  # print p-value for the test
  cat(paste0('The p-value of Chisq-test for testing the shape between X_t^2 and the mixture Z_t and Z_t^2 is ', round(p_value,4)))
  
  
  #####################
  # Count the portion of m^(1), m^(2)
  tree_pairwise <- lapply(tree_ranks_list, function(lst){
    lst[[1]]
  }) %>% 
    Reduce(f=rbind) %>% 
    pmin(5) %>% 
    as_tibble() %>% 
    `colnames<-`(c('x1', 'x2')) %>% 
    group_by(x1, x2) %>%
    summarise(n=n()) %>% 
    ungroup() %>% 
    mutate(freq = n/sum(n), type='tree')
  
  simulated_pairwise <- lapply(tree_ranks_list, function(lst){
    lst[[2]]
  }) %>% 
    Reduce(f=rbind) %>% 
    as.matrix() %>% 
    pmin(5) %>% 
    as_tibble() %>% 
    `colnames<-`(c('x1', 'x2')) %>% 
    group_by(x1, x2) %>%
    summarise(n=n()) %>% 
    ungroup() %>% 
    mutate(freq = n/sum(n), type='simulation')
  
  # Generate comparison bubble plot of the portions
  difference_plot <- rbind(tree_pairwise, simulated_pairwise) %>% 
    group_by(x1, x2) %>% 
    mutate(diff = freq[1]-freq[2]) %>% 
    filter(type=='tree') %>% 
    mutate(type= ifelse(diff>=0, 'tree', 'simulation')) %>% 
    ggplot(aes(x=x1, y=x2, size = abs(diff), colour=type)) + 
    scale_size(range = c(3, 10), name="Portion") +
    scale_colour_manual(values=c("#69b3a2", "#404080"), labels=lapply(c('$Simulation$', 'Tree Ring'), TeX), name="Type") +
    geom_point(alpha=0.3) +
    xlab(TeX('$m^{(1)}$')) +
    ylab(TeX('$m^{(2)}$')) +
    ggtitle(TeX('Pairwise Hermite Rank Portions Comparison')) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Confidence interval plot of P(m^(1)=1), P(m^(2)=1), P(m^(1)=1, m^(2)=1) for tree ring 
  se <- function(phat, n) {
    p_tilde <- (2 + n * phat) / (n + 4)
    return(sqrt(p_tilde * (1 - p_tilde) / n))
  }
  
  rank_comparison_summary <- lapply(c(1, 2), function(i) {
    tab <- lapply(tree_ranks_list, function(lst) {
      lst[[i]]
    }) %>%
      Reduce(f = rbind) %>%
      as_tibble() %>%
      `colnames<-`(c('x1', 'x2')) %>%
      group_by(x1, x2) %>%
      summarise(n12 = n()) %>%
      ungroup() %>%
      group_by(x1) %>%
      mutate(n1 = sum(n12)) %>%
      group_by(x2) %>%
      mutate(n2 = sum(n12)) %>%
      ungroup() %>%
      mutate(
        n = sum(n12),
        p12 = n12 / n,
        p1 = n1 / n,
        p2 = n2 / n,
        pcompare = sum(ifelse(x1<x2, n12, 0))/n[1]
      ) %>%
      filter(x1 == 1, x2 == 1) %>%
      select(pcompare, p1, p2, n)
  }) %>%
    Reduce(f = rbind) %>%
    mutate(Type = c('Tree Ring', 'Simulation')) %>%
    pivot_longer(!c('Type', 'n'), names_to = 'category', values_to = 'p') %>%
    mutate(se = ifelse(Type == 'Tree Ring', se(p, n), NA))
  
  pd <- position_dodge(0.2) # move them .05 to the left and right
  
  ci_plot <- rank_comparison_summary %>%
    ggplot(aes(
      x = factor(category),
      y = p,
      colour = Type
    )) +
    geom_errorbar(aes(ymin = p - 1.96 * se, ymax = p + 1.96 * se),
                  width = .2,
                  position = pd) +
    geom_point(position = pd) +
    scale_colour_manual(
      values = c("#69b3a2", "#404080"),
      labels = lapply(c('$Simulation$', 'Tree Ring'), TeX),
      name = "Type"
    ) +
    scale_x_discrete(
      limits = c('p1', 'p2', 'pcompare'),
      labels = c(
        'p1' = parse(text = TeX('$m^{(1)}=1$')),
        'p2' = parse(text = TeX('$m^{(2)}=1$')),
        'pcompare' = parse(text = TeX('$m^{(1)}\\, <\\,  m^{(2)}$'))
      )
    ) +
    ylim(c(0.35, 0.65)) +
    xlab(TeX('$Category$')) +
    ylab(TeX('$\\widehat{p}$')) +
    ggtitle(TeX(
      '$Point\\, Estimation\\,  and\\,  Confidence\\,  Interval\\,  of\\,  Proportions$'
    )) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # p <- grid_arrange_shared_legend(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], ncol = 2, nrow = 2, position='right')
  # plots <- grid.arrange(mix_comparison_plot, ci_plot,  ncol = 2)
  ggsave("plots/Squared_tree_ring_plot.png", plot = squared_tree_ring_plot,  width = 11/2, height = 7/2)
  ggsave("plots/Confidence_interval_plot.png", plot = ci_plot,  width = 11/2, height = 7/2)
  
}











# # Bubble plot for the portion of tree ring data for (m^(1), m^(2))
# portion_plot <- lapply(tree_ranks_list, function(lst){
#   lst[[1]]
# }) %>%
#   Reduce(f=rbind) %>%
#   pmin(5) %>%
#   as_tibble() %>%
#   `colnames<-`(c('x1', 'x2')) %>%
#   group_by(x1, x2) %>%
#   summarise(n=n()) %>%
#   ungroup() %>%
#   mutate(Portion = n/sum(n)) %>%
#   ggplot(aes(x=factor(x1), y=factor(x2), size = Portion)) +
#   geom_point(alpha=0.7) +
#   scale_x_discrete(
#     labels = c('1'='1', '2'='2', '3'='3', '4'='4','5'='4+'
#     )
#   ) +
#   scale_y_discrete(
#     labels = c('1'='1', '2'='2', '3'='3', '4'='4','5'='4+'
#     )
#   ) +
#   scale_size(range = c(3, 10), name="Portion") +
#   xlab(TeX('$m^{(1)}$')) +
#   ylab(TeX('$m^{(2)}$')) +
#   ggtitle(TeX('Pairwise Hermite Rank Portions of Tree Rings')) +
#   theme(plot.title = element_text(hjust = 0.5))


# lapply(tree_ranks_list, function(lst){
#   lst[[2]]
# }) %>% 
#   Reduce(f=rbind) %>% 
#   as.matrix() %>% 
#   pmin(5) %>% 
#   as_tibble() %>% 
#   `colnames<-`(c('x1', 'x2')) %>% 
#   group_by(x1, x2) %>%
#   summarise(n=n()) %>% 
#   ungroup() %>% 
#   mutate(freq = n/sum(n)) %>% 
#   ggplot(aes(x=x1, y=x2, size = freq)) +
#   geom_point(alpha=0.7)





# 
# 
# 
# table(data.frame(x=seq(1,3),y=seq(1,3))) 
# # simulated samples
# sapply(c(1,2), function(i){
#   tab <- lapply(tree_ranks_list, function(lst){
#     lst[[2]]
#   }) %>% 
#     Reduce(f=rbind) %>% 
#     as.matrix()
#   n <- dim(tab)[1]
#   tab_sum <- tab[,i] %>% 
#     table()
#   
#   p_hat <- tab_sum[1]/n
#   p_tilde <- (2+n*p_hat)/(n+4)
#   ci <- c(p_hat-1.96*sqrt(p_tilde*(1-p_tilde)/n), p_hat+1.96*sqrt(p_tilde*(1-p_tilde)/n))
#   
#   return(ci)
# }) %>% 
#   `colnames<-`(c('m1', 'm2'))
# # sapply(c(1,2), function(i){
# #   tab <- lapply(tree_ranks_list, function(lst){
# #     lst[[2]]
# #   }) %>% 
# #     Reduce(f=rbind) 
# #   nrow <- dim(tab)[1]
# #   tab_sum <- tab[,i] %>% 
# #     table()
# #   return(tab_sum[1]/nrow)
# # })
# 

