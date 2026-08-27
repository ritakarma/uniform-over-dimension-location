# Source required code files
source("code/utils.R")

# Load required packages
library("ggplot2")
library("patchwork")
library("colorspace")


# reading simulation results (output of the function 'run_simulation' 
# from the file "run_simulation.R")

r = readRDS("results/sim_results.rds")



# specify the valued of data dimension for plot
p_values = c(100, 50, 25, 5)

for (p in p_values) {
  # method names
  if (p == 5) {
    methods = c("KCDG2026^1", "KCDG2026^2", "sKCDG2026^1", "sKCDG2026^2", "ZGZC2020",
      "RS2016", "GBRSS2012[RBF]", "SKR2022[RBF]", "GBRSS2012[lin]", "SKR2022[lin]",
      "HT2", "CM1997"   
    )
  }
  else {
    methods = c("KCDG2026^1", "KCDG2026^2", "sKCDG2026^1", "sKCDG2026^2", "ZGZC2020",
       "RS2016", "GBRSS2012[RBF]", "SKR2022[RBF]", "GBRSS2012[lin]", "SKR2022[lin]",
       "BS1996", "CLX2014", "CQ2010", "SD2008", "CLZ2014"
    )
  }
  
  # define line types
  line_types = setNames(rep("solid", length(methods)), methods)
  line_types["KCDG2026^1"] <- "dashed"
  line_types["KCDG2026^2"] <- "dashed"
  line_types["sKCDG2026^1"] <- "dashed"
  line_types["sKCDG2026^2"] <- "dashed"
  
  # define line colors
  method_colors = setNames(
    qualitative_hcl(length(methods), h = c(0, 360), c = 100, l = 65),
    methods
  )
  
  # define legend labels
  
  legend_labels = parse(text = methods)
  
  
  # list for sub-plots 
  
  plots = list()
  
  # create the subplots
  
  for (i in 1:9) {
    # specify values of mu
    mu_values = 0:6
    if (p == 5) {
      mu_values = c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2)
    }
    
    plot_data = empirical_power(p = p, model = i, 
                                mu_values = mu_values, results = r)[ ,methods]
    plot_data = as.data.frame(as.table(plot_data))
    names(plot_data) = c("delta", "Method", "Value")
    plot_data$delta = as.numeric(sub("delta=", "", plot_data$delta))
    
    plots[[i]] <- ggplot(
      plot_data,
      aes(
        x = delta,
        y = Value,
        color = Method,
        linetype = Method,
        group = Method
      )
    ) +
      geom_line(linewidth = 0.5) +
      geom_point(size = 0.5) +
      
      scale_color_manual(
        values = method_colors,
        breaks = methods,
        labels = legend_labels
      ) +
      scale_linetype_manual(
        values = line_types,
        breaks = methods,
        guide = "none"
      ) +
      
      guides(
        color = guide_legend(
          override.aes = list(
            linetype = line_types[methods],
            size = 1.1
          )
        )
      ) +
      
      labs(
        title = paste("Model", new_label(paste0("Ex ", i))),
        x = expression(delta),
        y = "Power",
        color = NULL,
        linetype = NULL
      ) +
      
      theme_minimal() +
      
      theme(
        legend.position = "bottom",
        plot.title = element_text(
          hjust = 0.5,
          size = 10
        )
      )
  }
  
  # combine the subplots into a single plot
  
  final_plot <-
    wrap_plots(plots, ncol = 3) +
    plot_layout(guides = "collect") &
    theme(
      legend.position = "bottom",
      legend.margin = margin(r = 40),
      legend.key.width = unit(1.25, "cm"),
      legend.text = element_text(size = 8)
    )
  
  
  #print(final_plot)
  
  # create output file name
  file_name = paste0("figures/Rplot_p=", p, ".pdf")  
  
  # save the final plot
  ggsave(
    filename = file_name,
    plot = final_plot,
    width = 6.5,
    height = 8.25,
    units = "in"
  )
}


