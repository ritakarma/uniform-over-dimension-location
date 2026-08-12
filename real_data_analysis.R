# Colon tissue data is obtained from the package "plsgenomics".
# The dataset is also available at: 
# http://genomics-pubs.princeton.edu/oncology/affydata/index.html

library("plsgenomics")


# Source required code files

source("code/methods.R")

# Load required packages
library("HDNRA")
library("highmean")
library("maotai")
library("energy")

library("ggplot2")
library("patchwork")

## helper function to run all high dimensional tests

all_tests <- function(dataX, dataY, nsim = 10000, perm = 500) {
  # Function to perform all high dimensional tests
  #
  # Input: 
  #   dataX: n1 * p matrix -rows are X observations  
  #   dataY: n2 * p matrix -rows are Y observations
  #   nsim: number of simulations for cut-off calculation
  #   vec_beta:  vector specifying beta parameter values for second estimator
  #   perm: the number of permutation resamples to be used for kernel mmd and 
  #         energy tests, default value 500
  # Output: 
  #   a vector of p-values of different tests
  
  ## Performing our tests
  vec1 = kcdg_test(dataX = dataX, dataY = dataY, h = h_diff, 
                      nsim = nsim, 
                      vec_beta = 0.25)
  vec2 = kcdg_test(dataX = dataX, dataY = dataY, h = h_spatial, 
                      nsim = nsim, 
                      vec_beta = 0.25)
  vec = c(vec1, vec2)  
  
  names(vec) = c("KCDG2025^1", "KCDG2025^2", "sKCDG2025^1", "sKCDG2025^2")
  
  ## Performing other tests
  
  # Zhang et al 2020 Jasa (ZGZC2020)
  pval = HDNRA::ZGZC2020.TS.2cNRT(dataX, dataY)$p.value 
  vec = c(vec, ZGZC2020 = as.numeric(pval))  
  
  n1 = nrow(dataX)
  n2 = nrow(dataY)
  lab = c(rep(1,n1), rep(2,n2))  # group labels
  D = dist(rbind(dataX, dataY))  # distance matrix of pooled sample
  
  # Energy distance based test
  pval = energy::eqdist.etest(D, sizes = c(n1,n2), R = perm)$p.val 
  vec = c(vec, "RS2016" = as.numeric(pval))
  
  D = as.matrix(D)
  bw = median_bandwidth(D)  # median heuristic bandwidth for gaussian rbf kernel
  D = exp(-(D^2 / (2 * bw * bw)))   # gaussian rbf kernel
  
  # MMD based test with rbf kernel
  pval = maotai::mmd2test(K = D, 
                          label = lab, method = "u", mc.iter = perm)$p.value 
  vec = c(vec, "GBRSS2012[RBF]" = as.numeric(pval))
  
  # cross-MMD test with rbf kernel
  pval = crossMMD2sample(D, n1, n2)  
  vec = c(vec, "SKR2022[RBF]" = as.numeric(pval))
  
  
  D = tcrossprod(rbind(dataX, dataY)) # linear kernel
  
  # MMD based test with linear kernel
  pval = maotai::mmd2test(K = D, 
                          label = lab, 
                          method = "u", 
                          mc.iter = perm)$p.value 
  vec = c(vec, "GBRSS2012[lin]" = as.numeric(pval))
  
  # cross-MMD test with linear kernel
  pval = crossMMD2sample(D, n1, n2)  
  vec = c(vec, "SKR2022[lin]" = as.numeric(pval))
  
  # Bai and Saranadasa 1996 (BS1996)
  pval = highmean::apval_Bai1996(dataX, dataY)$pval 
  vec = c(vec, BS1996 = as.numeric(pval))
  
  # Cai et al 2014 (CLX2014)
  pval = highmean::apval_Cai2014(dataX, dataY)$pval 
  vec = c(vec, CLX2014 = as.numeric(pval))
  
  # Chen and Qin 2010 (CQ2010)
  pval = highmean::apval_Chen2010(dataX, dataY)$pval 
  vec = c(vec, CQ2010 = as.numeric(pval))
  
  # Chen et al 2014 (CLZ2014)
  pval = highmean::apval_Chen2014(dataX, dataY)$pval 
  vec = c(vec, CLZ2014 = as.numeric(pval))
  
  # Srivastava and Du 2008 (SD2008)
  pval = highmean::apval_Sri2008(dataX, dataY)$pval 
  vec = c(vec, SD2008 = as.numeric(pval))
  
  return(vec)
}


##############################################################################
##############################################################################


## Set seed for reporducibility
set.seed(100001)

# Loading Colon tissue data
data(Colon)
dataX = Colon$X[Colon$Y == 2, ]
dataY = Colon$X[Colon$Y == 1, ]
rm(Colon)

# perform tests on full dataset
result_data = all_tests(dataX = dataX, dataY = dataY)

# creating 50 data matrices by considering blocks of 40 consecutive 
# columns of the original data matrix
# i-th row of result_data_block contains the p-values of the tests 
# applied on the i-th block
result_data_block = c()
for (k in 1:50) {
  result_data_block = rbind(result_data_block, 
                            all_tests(dataX[ , (40 * (k-1) + 1):(40 * k)], 
                                      dataY[ , (40 * (k-1) + 1):(40 * k)])
                            )
}

result_data_avg = colSums(result_data_block) / nrow(result_data_block)


# layout of histogram subplots (number of histograms in each rows: 5, 5, 4) 
tests <- c(
  "ZGZC2020", "BS1996", "CLX2014", "CQ2010", "SD2008",
  "RS2016", "GBRSS2012[RBF]", "GBRSS2012[lin]", "SKR2022[RBF]", "SKR2022[lin]",
  "KCDG2025^1", "KCDG2025^2", "sKCDG2025^1", "sKCDG2025^2"
)

print("p-value of different tests for the whole dataset")
print(result_data[tests])

print("Average p-value over the data blocks for different tests")
print(result_data_avg[tests])


# Create a list of ggplot histograms
plots = list()
for (i in seq_along(tests)) {
  index = match(tests[i], colnames(result_data_block))
  # creating the i-th plot
  plots[[i]] <- ggplot(data.frame(p = result_data_block[, index]), aes(x = p)) +
    geom_histogram(
      breaks = seq(0, 1, by = 0.02),
      fill = "lightgray",
      color = "black"
    ) +
    
    scale_x_continuous(
      breaks = c(0, 0.5, 1)
      #labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0")
    )+
    
    labs(
      x = parse(text = tests[i]),
      y = "Frequency"
    ) +
    
    theme_minimal(base_size = 3) +
    
    theme(
      plot.margin = margin(1.5, 1.5, 1.5, 1.5),
      axis.title.x = element_text(size = 9, margin = margin(t = 2.5)),
      axis.title.y = element_text(size = 9, margin = margin(r = 2.5)),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

# Arrange in 3 rows: 5 on first row, 5 on second row, 4 on third row
row1 <- wrap_plots(plots[1:5], nrow = 1)
row2 <- wrap_plots(plots[6:10], nrow = 1)

row3 <- wrap_plots(
  plot_spacer(), plots[[11]], plots[[12]], plots[[13]], plots[[14]], plot_spacer(),
  nrow = 1
) + plot_layout(widths = c(0.5, 1, 1, 1, 1, 0.5))

final_plot <- row1 / row2 / row3


# save final plot

ggsave(
  filename = "figures/Rplot_histogram.pdf",
  plot = final_plot,
  width = 6.5,
  height = 4.5,
  units = "in"
)


  