# Source required code files
source("code/utils.R")

library(knitr)


# reading simulation results (output of the function 'run_simulation' 
# from the file "run_simulation.R")

r = readRDS("results/sim_results.rds")


p_values = c(100, 50, 25, 5)

for (p in p_values) {
  ## preparing data frame to print as a table
  m = r[[paste0("p=", p, "_mu=0")]]$out
  rownames(m) = new_label_tex(new_label(rownames(m)))
  colnames(m) = new_label(colnames(m))
  
  # changing order of rows for final table
  if (p == 5) {
    m <- m[c(1:5, 11, 12, 9, 10, 7, 8, 6), , drop = FALSE]
  }
  else {
    m <- m[c(1:5, 11:15, 9, 10, 7, 8, 6), , drop = FALSE] 
  }
  
  m = as.data.frame(m, check.names = FALSE)
  m = cbind(Model = rownames(m), m)
  rownames(m) = NULL
  
  ## printing the data frame in console
  print(paste0("Estimated sizes at nominal level 5% for different tests for p = ", p))
  print(m)
  
  # creating latex code for the table
  tex <- kable(
    m,
    format = "latex",
    col.names = colnames(m),
    align = rep("c", ncol(m)),
    booktabs = TRUE,
    escape = FALSE
  )
  
  tex <- gsub("\\\\toprule|\\\\midrule|\\\\bottomrule", "\\\\hline", tex)
  tex <- gsub("\\\\addlinespace", "", tex)
  
  writeLines(tex, paste0("tables/Rtable_p=", p, ".tex"))
}



