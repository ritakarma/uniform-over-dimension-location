#### Requires: "HDNRA", "mvtnorm", "highmean", "DescTools", "maotai", "energy"
#### Requires: "code/methods.R"

model_sim <- function(index, n1, n2, p, mu) {
  # Generates data for simulation study
  #
  # Input: 
  #  index = model index (between 1 to 9)
  #  n1 = number of X observations
  #  n2 = number of Y observations
  #  p = dimension
  #  mu = parameter for location shift
  #       [amount of location shift is c * mu * h_vec where 
  #        h_vec = (1:p)/sqrt(sum((1:p)^2)) and c = 10 for models 3, 6, and 9
  #        and c = 1 for the remaining models]
  # Output: 
  #  dataX: n1 * p matrix 
  #  dataY: n2 * p matrix
  
  h_vec = (1:p)/sqrt(sum((1:p)^2)) ## location shift of Y from X 
  ## is proportional to mu*h_vec
  
  ## Models
  
  ##Model 1
  if (index == 1) {    
    ## X: N(0, Sigma), Y: N(0, Sigma), Sigma: equicorrelation matrix with 
    ## non-diagonal elements 0.5
    
    Sigma = matrix(0.5, nrow = p, ncol = p)
    diag(Sigma) = rep(1, p)
    dataX = mvtnorm::rmvnorm(n1, mean = rep(0, p), sigma = Sigma)
    dataY = mvtnorm::rmvnorm(n2, mean = mu*h_vec, sigma = Sigma)
    return(list(dataX = dataX, dataY = dataY))
  }
  
  ## Model 2
  if (index == 2) {
    ## X: t_4(0, Sigma), Y: t_4(0, Sigma), Sigma: equicorrelation matrix with 
    ## non-diagonal elements 0.5
    
    Sigma = matrix(0.5, nrow = p, ncol = p)
    diag(Sigma) = rep(1, p)
    dataX = mvtnorm::rmvt(n1, sigma = Sigma, df = 4, 
                          delta = rep(0, p), type = "shifted")
    dataY = mvtnorm::rmvt(n2, sigma = Sigma, df = 4, 
                          delta = mu*h_vec, type = "shifted")
    return(list(dataX = dataX, dataY = dataY))
  }
  
  ## Model 3
  if (index == 3) {
    ## X: t_1(0, Sigma), Y: t_1(0, Sigma), Sigma: equicorrelation matrix with
    ## non-diagonal elements 0.5
    
    Sigma = matrix(0.5, nrow = p, ncol = p)
    diag(Sigma) = rep(1, p)
    dataX = mvtnorm::rmvt(n1, sigma = Sigma, df = 1, 
                          delta = rep(0, p), type = "shifted")
    dataY = mvtnorm::rmvt(n2, sigma = Sigma, df = 1, 
                          delta = 10*mu*h_vec, type = "shifted")
    return(list(dataX = dataX, dataY = dataY))
  }
  
  ##Model 4
  if (index == 4) {
    ## X: N(0, Sigma), Y: N(0, Sigma), Sigma: identity
    
    Sigma = diag(1, nrow = p)
    dataX = mvtnorm::rmvnorm(n1, mean = rep(0, p), sigma = Sigma)
    dataY = mvtnorm::rmvnorm(n2, mean = mu*h_vec, sigma = Sigma)
    return(list(dataX = dataX, dataY = dataY))
  }
  
  ## Model 5
  if (index == 5) {
    ## X: t_4(0, Sigma), Y: t_4(0, Sigma), Sigma: identity
    
    Sigma = diag(1, nrow = p)
    dataX = mvtnorm::rmvt(n1, sigma = Sigma, df = 4, 
                          delta = rep(0, p), type = "shifted")
    dataY = mvtnorm::rmvt(n2, sigma = Sigma, df = 4, 
                          delta = mu*h_vec, type = "shifted")
    return(list(dataX = dataX, dataY = dataY))
  }
  
  ## Model 6
  if (index == 6) {
    ## X: t_1(0, Sigma), Y: t_1(0, Sigma), Sigma: identity
    
    Sigma = diag(1, nrow = p)
    dataX = mvtnorm::rmvt(n1, sigma = Sigma, df = 1, 
                          delta = rep(0, p), type = "shifted")
    dataY = mvtnorm::rmvt(n2, sigma = Sigma, df = 1, 
                          delta = 10*mu*h_vec, type = "shifted")
    return(list(dataX = dataX, dataY = dataY))
  }
  
  ##Model 7
  if (index == 7) {
    ## X: N(0, Sigma), Y: N(0, Sigma), Sigma: (i,j)-th entry (0.7)^|i-j|
    
    DD = abs(outer(1:p, 1:p, "-"))
    Sigma = (0.75)^DD
    dataX = mvtnorm::rmvnorm(n1, mean = rep(0, p), sigma = Sigma)
    dataY = mvtnorm::rmvnorm(n2, mean = mu*h_vec, sigma = Sigma)
    return(list(dataX = dataX, dataY = dataY))
  }
  
  ## Model 8
  if (index == 8) {
    ## X: t_4(0, Sigma), Y: t_4(0, Sigma), Sigma: (i,j)-th entry (0.7)^|i-j|
    
    DD = abs(outer(1:p, 1:p, "-"))
    Sigma = (0.75)^DD
    dataX = mvtnorm::rmvt(n1, sigma = Sigma, df = 4, 
                          delta = rep(0, p), type = "shifted")
    dataY = mvtnorm::rmvt(n2, sigma = Sigma, df = 4, 
                          delta = mu*h_vec, type = "shifted")
    return(list(dataX = dataX, dataY = dataY))
  }
  
  ## Model 9
  if (index == 9) {
    ## X: t_1(0, Sigma), Y: t_1(0, Sigma), Sigma: (i,j)-th entry (0.7)^|i-j|
    
    DD = abs(outer(1:p, 1:p, "-"))
    Sigma = (0.75)^DD
    dataX = mvtnorm::rmvt(n1, sigma = Sigma, df = 1, 
                          delta = rep(0, p), type = "shifted")
    dataY = mvtnorm::rmvt(n2, sigma = Sigma, df = 1, 
                          delta = 10*mu*h_vec, type = "shifted")
    return(list(dataX = dataX, dataY = dataY))
  }
  
}



simulate_tests <- function(n1, n2, p, mu, models, Rep) {
  # Runs experiments using simulated data
  #
  # Input: 
  #  n1 = number of X observations
  #  n2 = number of Y observations
  #  p = dimension
  #  mu = parameter for location shift
  #  models = indices of models generated using model_sim
  #  Rep = number of replications to calculate the proportion of rejections
  # Output: 
  #  a matrix with entries equal to the proportion of rejections and 
  #  columns and rows corresponding to the models and tests respectively
  
  
  ## specify fixed paramaters below
  
  ## specify level of the test
  alpha = 0.05
  
  ## specify the number of permutation resamples to be used 
  ## for kernel mmd and energy tests
  perm = 500
  
  ## tests corresponding to tapering estimators for different beta's
  Beta = 0.25 #c(0.1,0.2,0.3,0.4,0.5)  
    
  if (length(Beta) > 0) {
    tests = c("KCDG^1",
              paste0("KCDG^2_", Beta), "sKCDG^1", paste0("sKCDG^2_", Beta))
  }
  else {
    tests = c("KCDG^1", "sKCDG^1")
  }
  
  tests = c(tests, "ZGZC2020", "Energy",
            "MMD_rbf", "crossMMD_rbf", "MMD_linear", "crossMMD_linear")
  
  ## Other tests
  if (p <= 10) {
    tests = c(tests, "T^2", "CM1997")
  }
  else {
    tests = c(tests, "BS1996", "CLX2014", "CQ2010", "CLZ2014", "SD2008") 
  }
  
  results = matrix(0, nrow = length(tests), ncol = length(models),
                   dimnames = list(tests, paste0("Ex ", as.character(models))))
  
  for (id in 1:length(models)) {
    vec_final = numeric(length(tests))  # vector to store total 
                                        # number of rejections
    
    for (i  in 1:Rep) {
      res = model_sim(models[id], n1, n2, p, mu) ## Generating dataX and dataY
      dataX = res$dataX
      dataY = res$dataY
      
      ## Performing our test with kernel h_diff, h_spatial and specified beta's

      vec1 = kcdg_test(dataX = dataX, dataY = dataY,
                          h = h_diff, nsim = 10000, vec_beta = Beta)
      vec2 = kcdg_test(dataX = dataX, dataY = dataY,
                          h = h_spatial, nsim = 10000, vec_beta = Beta)
      vec = c(vec1, vec2)
      
      ## Performing other tests
      
      # Zhang et al 2020 Jasa (ZGZC2020)
      pval = HDNRA::ZGZC2020.TS.2cNRT(dataX, dataY)$p.value 
      vec = c(vec, as.numeric(pval))   
      
      lab = c(rep(1,n1), rep(2,n2))  # group labels
      D = dist(rbind(dataX, dataY))  # distance matrix of pooled sample
      
      # Energy distance based test
      pval = energy::eqdist.etest(D, sizes = c(n1,n2), R = perm)$p.val 
      vec = c(vec, as.numeric(pval))
      
      D = as.matrix(D)
      bw = median_bandwidth(D)          # median heuristic bandwidth for 
                                        # gaussian rbf kernel
      D = exp(-(D^2 / (2 * bw * bw)))   # gaussian rbf kernel
      
      # MMD based test with rbf kernel
      pval = maotai::mmd2test(K = D, 
                              label = lab,
                              method = "u",
                              mc.iter = perm)$p.value 
      vec = c(vec, as.numeric(pval))
      
      
      # cross-MMD test with rbf kernel
      pval = crossMMD2sample(D, n1, n2)  
      vec = c(vec, as.numeric(pval))
      
      
      D = tcrossprod(rbind(dataX, dataY)) # linear kernel
      
      
      # MMD based test with linear kernel
      pval = maotai::mmd2test(K = D, 
                              label = lab, 
                              method = "u",
                              mc.iter = perm)$p.value 
      vec = c(vec, as.numeric(pval))
      
      
      # cross-MMD test with linear kernel
      pval = crossMMD2sample(D, n1, n2)  
      vec = c(vec, as.numeric(pval))
      
      if (p <= 10) {
        # Hotelling's T^2 test
        pval = DescTools::HotellingsT2Test(x = dataX, y = dataY)$p.value 
        vec = c(vec, as.numeric(pval))
        
        # Choi and Marden (1997)
        pval = choi_marden_2sample(dataX, dataY, h = h_spatial) 
        vec = c(vec, as.numeric(pval))
      }
      else {
        # Bai and Saranadasa 1996 (BS1996)
        pval = highmean::apval_Bai1996(dataX, dataY)$pval 
        vec = c(vec, as.numeric(pval))
        
        # Cai et al 2014 (CLX2014)
        pval = highmean::apval_Cai2014(dataX, dataY)$pval 
        vec = c(vec, as.numeric(pval))
        
        # Chen and Qin 2010 (CQ2010)
        pval = highmean::apval_Chen2010(dataX, dataY)$pval 
        vec = c(vec, as.numeric(pval))
        
        # Chen et al 2014 (CLZ2014)
        pval = highmean::apval_Chen2014(dataX, dataY)$pval 
        vec = c(vec, as.numeric(pval))
        
        # Srivastava and Du 2008 (SD2008)
        pval = highmean::apval_Sri2008(dataX, dataY)$pval 
        vec = c(vec, as.numeric(pval)) 
      }
      
      # reject null when p-value is less than alpha
      vec = as.numeric(vec <= alpha)
      
      vec_final = vec_final + vec
    }
    results[,id] = vec_final / Rep  # proportion of rejections
  }
  
  return(results)
}