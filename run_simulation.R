# Source required code files
source("code/methods.R")
source("code/simulation.R")

# Load required packages
library("HDNRA")
library("mvtnorm")
library("highmean")
library("DescTools")
library("maotai")
library("energy")

library("future")
library("future.apply")

# helper function to run the simulation study
run_simulation <- function(p_values, mu_values, Rep) {
  # Runs the complete simulation study using the given p and mu values
  #
  # Input:
  #   p_values: vector of values of data dimension (must be positive integers)
  #   mu_values: vector of possible values of mu (parameters for location shift)
  #              see 'mode_sim' function from the file 'simulation.R' for details
  #   Rep: number of replications to calculate the proportion of rejections
  # Returns:
  #    Named list. Each element corresponds to a unique (p, mu). Names have the 
  #    form "p={p}_mu={mu}", e.g. "p=100_mu=3".
  #    parameter combination and is itself a named list containing:
  #      p: Value of p.
  #      mu: Value of mu.
  #      seed: Deterministic seed used to generate the output.
  #      out: Named matrix with row names corresponding to the testing method
  #      (e.g KCDG^1, ZGZC2020 etc.) and column names corresponding to the 
  #      experiment number (e.g "EX 1"), and the entries are the proportion of 
  #      times the the null hypothesis was rejected (i.e. the empirical 
  #      probability of type 1 error if mu=0 and the empirical power when mu>0)
  
  ## Fix parameter values
  
  # Specify sample size
  n1 = 40
  n2 = 50
  
  ## Fix data generating models
  models = 1:9
  
  grid = expand.grid(mu = unique(mu_values), p = unique(p_values))
  
  results = future_lapply(
    X = seq_len(nrow(grid)),
    FUN = function(i){
      # Current parameter combination
      p = grid$p[i]
      mu = grid$mu[i]
      
      # Deterministic seed for this combination
      
      seed = 100 * p + mu
      set.seed(seed)
      
      # Run simulation
      res = simulate_tests(n1, n2, p, mu, models, Rep)
      
      return(list(p = p, mu = mu, seed = seed, out = res))
    },
    future.seed = NULL
  )
  
  names(results) = paste0("p=", grid$p, "_mu=", grid$mu)
  
  return(results)
}



#### perform the complete simulation study

start = Sys.time()

## number of replications to calculate the proportion of rejections

Rep = 10000

## Set-up parallel processing

num_workers = max(1, availableCores() - 1)
old_plan = plan(multisession, workers = num_workers)

## specify values of p and mu (high/medium dimensional setup)
p_values = c(100, 50, 25)
mu_values = 0:6 

## run simulation study
results = run_simulation(p_values = p_values, mu_values = mu_values, Rep = Rep)

## specify values of p and mu (low dimensional setup)
p_values = 5
mu_values = c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2)

## run simulation study and combine results
results = c(results,
          run_simulation(p_values = p_values, mu_values = mu_values, Rep = Rep)
)

# restore initial plan
plan(old_plan)

### save output in a rds file

saveRDS(results, 
        file = "results/sim_results.rds")

end = Sys.time()

print((end - start))



