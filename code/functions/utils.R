new_label <- function(x) {
  # function to modify method / experiment  labels for final plot / table
  #
  # Input: 
  #   character vector with names / labels of methods / experiments
  # Output: 
  #   character vector with modified names
  
  map = c("KCDG^1" = "KCDG2026^1", "KCDG^2_0.25" = "KCDG2026^2",
          "sKCDG^1" = "sKCDG2026^1", "sKCDG^2_0.25" = "sKCDG2026^2", 
          "T^2" = "HT2",
          "MMD_rbf" = "GBRSS2012[RBF]", "MMD_linear" = "GBRSS2012[lin]",
          "crossMMD_rbf" = "SKR2022[RBF]", "crossMMD_linear" = "SKR2022[lin]",
          "Energy" = "RS2016",
          "Ex 1" = "1.i.", "Ex 2" = "2.i.", "Ex 3" = "3.i.",
          "Ex 4" = "1.ii.", "Ex 5" = "2.ii.", "Ex 6" = "3.ii.",
          "Ex 7" = "1.iii.", "Ex 8" = "2.iii.", "Ex 9" = "3.iii."
  )
  matched = x %in% names(map)
  x[matched] = map[x[matched]]
  return(x)
}



empirical_power <- function(p, model, mu_values, results) {
  # Takes the simulation outcome as input and returns a matrix containing 
  # the empirical power / size of the different tests under the specified
  # dimension, and experiment / model number and values of mu.
  #
  # Input: 
  #   p: value of dimension p
  #   model: experiment/model number
  #   mu_values: vector of values of mu
  #     In our case it is 
  #       - 0:6 for p = 100, 25, 50
  #       - c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2) for p = 5 
  #   results: result of the simulation study (output of the function
  #            'run_simulation' from the file "run_simulation.R")
  # Returns:
  #   A named matrix with rows indexed by values (e.g. "delta=0",
  #   "delta=0.2" etc.), 
  #   and columns indexed by model/expriment (e.g "Ex 1", "Ex 2" etc.),
  #   the entries are the proportion of times the the null hypothesis 
  #   was rejected (i.e. the empirical probability of type 1 error if mu=0 
  #   and the empirical power when mu > 0)
  
  rows = list()
  col_name = paste0("Ex ", model)
  row_names = c()
  
  # Location shift is delta * h_vec where h_vec = (1:p)/sqrt(sum((1:p)^2)) 
  # delta = c * mu where c = 10 for models 3, 6, and 9 and c = 1 for the 
  # remaining models (see the function 'model_sim' in "simulation.R")
  
  c = 1
  if (model %in% c(3, 6, 9)) {
    c = 10
  }
  
  for (i in seq_along(mu_values)) {
    m = results[[paste0("p=", p, "_mu=", mu_values[i])]]$out
    rows[[i]] = m[,col_name]
    if (!is.null(m)){
      row_names = c(row_names, paste0("delta=", c*mu_values[i]))
    }
  }
  m = do.call(rbind, rows)
  rownames(m) = row_names
  colnames(m) = new_label(colnames(m))
  return(m)
}


new_label_tex <- function(x) {
  # function to modify method / experiment  labels for final plot / table
  #
  # Input: 
  #   character vector with names / labels of methods / experiments
  # Output: 
  #   character vector with modified names
  
  map = c("KCDG2026^1" = "KCDG2026$^1$", "KCDG2026^2" = "KCDG2026$^2$", 
          "sKCDG2026^1" = "sKCDG2026$^1$", "sKCDG2026^2" = "sKCDG2026$^2$",
          "GBRSS2012[RBF]" = "GBRSS2012$_{\\text{RBF}}$", 
          "GBRSS2012[lin]" = "GBRSS2012$_{\\text{lin}}$",
          "SKR2022[RBF]" = "SKR2022$_{\\text{RBF}}$", 
          "SKR2022[lin]" = "SKR2022$_{\\text{lin}}$"
  )
  matched = x %in% names(map)
  x[matched] = map[x[matched]]
  return(x)
}

