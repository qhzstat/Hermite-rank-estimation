library(EQL)
library(latex2exp)

source("functions\\parallel.R")
source("functions\\LRD_Generator_new.R")

if(sys.nframe()==0){
  alpha <- 0.2
  n <- 5000
  nrep <- 3000
  set.seed(44)
  
  G_list <- list()
  G_list[[1]] <- function(x) hermite(x, n=1, prob = T) 
  G_list[[2]] <- function(x) hermite(x, n=2, prob = T) 
  G_list[[3]] <- function(x) hermite(x, n=3, prob = T) 
  G_list[[4]] <- function(x) hermite(x, n=4, prob = T) 
  
  hist_plots <- mclapply(seq(1,4), function(i){
    G <- G_list[[i]]
    xbars <- LRD_gen(n=n, alpha=alpha, nrep=nrep, option='FGaussian') %>% 
      sapply(function(zt) mean(G(zt)))
    xbars <- xbars/sd(xbars)
    hist_plot <- tibble(xbar = xbars) %>% 
      ggplot(aes(x=xbar)) +
      geom_histogram(aes(y=..density..), colour="black", fill="white") +
      geom_density(alpha=.2, fill="#69b3a2") +
      ylim(c(0, 0.9)) + 
      xlim(c(-8,8)) + 
      xlab(TeX('$\\bar{X}_n/\\sqrt{Var(\\bar{X}_n)}$')) + 
      ggtitle(TeX(paste0('$Histogram\\, (X_t=H_', i, '(Z_t))$'))) +
      theme(plot.title = element_text(hjust = 0.5))
    
    return(hist_plot)
  })
  
  hist_xbars <- grid.arrange(hist_plots[[1]], hist_plots[[2]], hist_plots[[3]], hist_plots[[4]], nrow = 2)
  
  ggsave("plots/Histogram_xbar.png", plot = hist_xbars,  width = 11, height = 7)
}


