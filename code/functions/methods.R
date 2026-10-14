h_diff <- function(dataX, dataY) {
  # Computes array with entried h(X_i, X_j), h being the difference kernel
  # for KCDG or sKCDG tests
  #
  # Input: 
  #   dataX: n1 * p matrix -rows are X observations  
  #   dataY: n2 * p matrix -rows are Y observations
  # Returns: 
  #   array of size n1 * n2 * p with (i, j, )-th entry X_i - Y_j

  n1 = nrow(dataX)
  n2 = nrow(dataY)
  p = ncol(dataX)
  X = array(dataX, dim = c(n1, p, n2))  # n1 * p * n2 array where each slice
                                        # ( , , k) is a copy of dataX
  
  X = aperm(X, c(1, 3, 2))              # swapped dimensions to 
                                        # make n1 * n2 * p array  
  
  Y = array(dataY, dim = c(n2, p, n1))  # n2 * p * n1 array where each 
                                        # slice ( , , k) is a copy of dataY
  
  Y = aperm(Y, c(3, 1, 2))              # swapped dimensions to make 
                                        # n1 * n2 * p array  
  
  return(X - Y)
}




h_spatial <- function(dataX, dataY, tol = 1e-8) {
  # Computes array with entried h(X_i, X_j), h being the spatial
  # kernel for KCDG or sKCDG tests
  #
  # Input: 
  #   dataX: n1 * p matrix -rows are X observations  
  #   dataY: n2 * p matrix -rows are Y observations
  #   tol: tolerance to check if a number is zero
  # Returns: 
  #   array of size n1 * n2 * p with (i, j, )-th entry 
  #   (X_i - Y_j) / ||X_i - Y_j|| if ||X_i - Y_j|| > 0 and 0 otherwise
  
  n1 = nrow(dataX)
  n2 = nrow(dataY)
  p = ncol(dataX)
  X = array(dataX, dim = c(n1, p, n2))    # n1 * p * n2 array where each slice
                                          # ( , , k) is a copy of dataX
  
  X = aperm(X, c(1, 3, 2))      # swapped dimensions to make n1 * n2 * p array  
  
  Y = array(dataY, dim = c(n2, p, n1))    # n2 * p * n1 array where each 
                                          # slice ( , , k) is a copy of dataY
  
  Y = aperm(Y, c(3, 1, 2))      # swapped dimensions to make n1 * n2 * p array  
  
  D = X - Y                 # n1 * n2 * p array with (i, j, )-th entry X_i - Y_j
  
  # swap columns and then sum over 1st dimension 
  weights = colSums(aperm((D * D), c(3, 1, 2)))  # n1 * n2 matrix with (i,j)-th
                                                 # entry ||X_i - Y_j||^2
  
  weights = sqrt(weights)
  id = (weights > tol)
  
  weights[id] = 1 / weights[id]     # (i,j)-th entry reciprocal of ||X_i - Y_j||
  weights[!id] = 0                  # (i,j)-th entry 0 if ||X_i - Y_j|| <= tol
  
  return(D * array(weights, dim = c(n1, n2, p)))
}





compute_T_Sigma <- function(dataX, dataY, h) {
  # Computes the test statistic for KCDG or sKCDG tests
  #
  # Input: 
  #   dataX: n1 * p matrix -rows are X observations  
  #   dataY: n2 * p matrix -rows are Y observations
  #   h: kernel function - should take inputs dataX, dataY and give output an 
  #      array of size n1 * n2 * p with (i, j, )-th entry h(X_i, Y_j)
  # Returns: 
  #   A list containing -
  #    T: the test statistic T_{n,p}
  #    Sigma: estimator of Sigma_p  
  
  n1 = nrow(dataX)
  n2 = nrow(dataY)
  p = ncol(dataX)
  
  #if (p != ncol(dataY)) {
  #  stop("Datasets have different dimensions")
  #}
  #if (n1 <= 1 || n2 <= 1 || p == 0) {
  #  stop("Must have rowsize of both datasets at least 2 and columnsize at least 1")
  #}
  
  n = n1 + n2
  H = h(dataX, dataY)     #dim: n1 * n2 * p
  
  if (!is.array(H) || !all(dim(H) == c(n1, n2, p))) {
    stop("Output of 'h' must be of dimension: 
         n1 * n2 * p where n1 = nrow(dataX), n2 = nrow(dataY),
         p = ncol(dataX) = ncol(dataY)")
  }
  
  
  # sums over 1st dimension of H
  S2 = colSums(H)                    # S2: (n2 * p matrix), j-th row is 
                                     # \sum_{i=1}^{n_2} h(X_i, Y_j)^t
  
  # swaps 1st and 2nd dimension and then sums over 1st dimension of H
  S1 = colSums(aperm(H, c(2,1,3)))   # S1: (n1 * p matrix), i-th row is 
                                     # \sum_{j=1}^{n_2} h(X_i, Y_j)^t
  
  
  S = colSums(S1)             # p * 1 matrix/vector: sum_i sum_j h(X_i, Y_j) 
  
  
  S = tcrossprod(S)      # p * p matrix: 
                         # (sum_i sum_j h(X_i, Y_j)) (sum_i sum_j h(X_i, Y_j))^t 
  
  S1 = crossprod(S1)     # p * p matrix: 
                         # sum_i (sum_j h(X_i, Y_j)) (sum_j h(X_i, Y_j))^t
  
  S2 = crossprod(S2)     # p * p matrix: 
                         # sum_j (sum_i h(X_i, Y_j)) (sum_i h(X_i, Y_j))^t
  
  
  sum_S_sq = sum(diag(S))     # scalar: ||sum_i sum_j h(X_i, Y_j)||^2
  
  sum_S1_sq = sum(diag(S1))   # scalar: sum_i ||sum_j h(X_i, Y_j)||^2
  
  sum_S2_sq = sum(diag(S2))   # scalar: sum_j ||sum_i h(X_i, Y_j)||^2
  
  sum_sq_h = sum(H * H)       # scalar: sum_i sum_j ||h(X_i, Y_j)||^2
  
  
  T_denom = n * n1 * n2
  
  T_val = ( sum_S_sq - sum_S1_sq - sum_S2_sq + sum_sq_h ) / T_denom
  
  Sigma = ( S1 + S2 ) / T_denom - S / (n1 * n1 * n2 * n2)
  
  return(list(T = T_val, Sigma = Sigma))
}



compute_Sigma_tap <- function(Sigma, n, beta) {
  # Computes tapered covariance matrix required for the implementation of
  # sKCDG test
  #
  # Input: 
  #   Sigma: p * p matrix Sigma_p
  #   n: (n1 + n2) - total number of observations
  #   beta: tuning parameter beta > 0
  # Returns: 
  #   tapered version of the input matrix (used to calculate \hat\Sigma_p^{(2)}) 
  
  p = nrow(Sigma)
  
  pow = 1 / (2 * beta + 2)
  k = floor(n^pow)                                     
  if (k > p) {
    k = p
  }                                 # k = min{p, floor(n^(1 / (2 * beta + 2)))}
  
  D = abs(outer(1:p, 1:p, "-"))     # D_(ij) = |i-j|
  
  W = matrix(0, nrow = p, ncol = p)
  
  W[D <= (k/2)] = 1                 # W_(ij) = 1 if |i-j| <= k/2
  
  middle_indices = (D > (k/2) & D < k)
  
  # W_(ij) = 2 * (1 - D_(ij) / k) if |i-j| > k/2 and |i-j| < k
  W[middle_indices] = 2 * (1 - D[middle_indices] / k)  
  
  return (Sigma * W)
}





kcdg_test <- function(dataX, 
                      dataY, 
                      h = h_spatial, 
                      nsim, 
                      estimators = c(1,1), 
                      vec_beta = 0.25, 
                      batchsize = 10000) {
  # Performs the KCDG or sKCDG test with plain and/or tapering estimator
  # and returns a vector of p-value(s)
  #
  # Input: 
  #   dataX: n1 * p matrix - rows are X observations  
  #   dataY: n2 * p matrix - rows are Y observations
  #
  #   h: kernel function [this function should take inputs dataX, dataY and 
  #                       give output an array of size n1 * n2 * p with
  #                       (i, j, )-th entry equal to h(X_i, Y_j)]
  #     - use h_diff for KCDG test
  #     - use h_spatial for sKCDG test
  #
  #   nsim: number of simulations for cut-off calculation
  #
  #   estimators: a vector of length 2 - if 1st entry is 1 plain estimator is 
  #               used and if 2nd entry is 1 tapering estimator is used, both
  #               estimators are used when the argument is c(1,1)
  #
  #   vec_beta:  vector specifying beta parameter values for second estimator
  #
  #   batchsize: size of batches for simulation
  #
  # Returns: 
  #   A vector of p-values with names given by: 
  #    "plain": corresponds to plain estimator
  #    "tap_beta={beta}": corresponds to tapering estimator with beta parameter 
  #                       value equal to {beta}
  
  if (all(estimators == 0)) {
    return(NULL)
  }
  
  n1 = nrow(dataX)
  n2 = nrow(dataY)
  p1 = ncol(dataX)
  p2 = ncol(dataY)
  
  ## sanity check
  if (n1 < 2) stop("Number of rows in first data matrix is < 2")
  if (n2 < 2) stop("Number of rows in second data matrix is < 2")
  if (p1 != p2) stop("Two data matrices have a different number of columns")

  
  res = compute_T_Sigma(dataX, dataY, h)
  T = res$T
  n = n1 + n2
  
  # vector of indices corresponding to the estimators
  if (estimators[1] == 1 & estimators[2] == 1 & length(vec_beta) > 0) {
    indices = 0:length(vec_beta)
    vec_name = c("plain", paste0("tap_beta=", as.character(vec_beta))) 
  }
  else if (estimators[1] == 1){
    indices = 0
    vec_name = "plain"
  }
  else {
    indices = 1:length(vec_beta)
    vec_name = paste0("tap_beta=", as.character(vec_beta))
  }
  
  # list of p-values
  results = numeric(length(indices))
  names(results) = vec_name
  
  for (i in 1:length(indices)) {
    if (indices[i] == 0) {
      Sigma = res$Sigma  # plain estimator
    }
    else {
      #tapering estimator
      Sigma = compute_Sigma_tap(res$Sigma, n, vec_beta[indices[i]])  
    }
    
    # vector of eigenvalues
    lambda <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
    lambda <- pmax(lambda, 0)  # clip small negative eigenvalues if exists
    p = length(lambda)
    
    # vector of simulated random variables
    Q = numeric(nsim)
    
    # number of batches
    nbatches = ceiling(nsim / batchsize)
    
    # loop over batches
    for (b in 1:nbatches) {
      idx = ((b - 1) * batchsize + 1) : (min(b * batchsize, nsim))
      
      # generate standard normals
      Z = matrix(rnorm(length(idx) * p), nrow = length(idx), ncol = p)
      
      # compute quadratic forms for this batch
      Q[idx] = (Z^2 - 1) %*% lambda
    }
    
    # p-value
    results[i] = mean(Q >= T) 
  }
  
  return(results)
}



choi_marden_2sample <- function(dataX, dataY, h = h_spatial) {
  # Perfoms two sample test by Choi and Marden (1997) and returns p-value
  # Paper Link: 
  #    https://www.tandfonline.com/doi/pdf/10.1080/01621459.1997.10473680
  #
  # Input: 
  #   dataX: n1 * p matrix -rows are X observations  
  #   dataY: n2 * p matrix -rows are Y observations
  #   h: kernel function (should take inputs dataX, dataY and give output an 
  #      array of size n1 * n2 * p with (i, j, )-th entry h(X_i, Y_j))
  #      -- h = h_diff: Hotelling's T^2 test with chi-square approximation
  #      -- h = h_spatial (default): Spatial Rank test by Choi and Marden (1997) 
  # Returns: 
  #   p-value of the test
  
  n1 = nrow(dataX)
  n2 = nrow(dataY)
  n = n1 + n2
  p = ncol(dataX)
  
  #### Calculating Sigma_hat
  
  # h(dataX, dataX) is n1 * n1 * p array with (i,j, )-th entry being h(X_i, X_j)
  # S is n1 * p matrix with j-th row - 
  # \sum_{i=1}^{n_1} h(X_i, X_j)^t / n1 =  -R_N1(X^(j))^t 
  # as defined in Choi and Marden (1997)
  S = colSums(h(dataX, dataX)) / n1
  
  # S^t S = \sum_{j=1}^{n_1} R_N1(X^(j)) R_N1(X^(j))^t
  Sigma_hat = crossprod(S) 
  
  # Same calculation with the second group
  S = colSums(h(dataY, dataY)) / n2
  Sigma_hat = Sigma_hat + crossprod(S)
  Sigma_hat = Sigma_hat / (n-2)
  
  #### Calculating test statistic KW
  # n1*n*\bar{R}^(1) 
  # = \sum_{i=1}^n1 n * R(X^(i)) 
  # = \sum_{i=1}^n1 [\sum_i' h(X_i, X_i') + \sum_j h(X_i, Y_j)]
  # = \sum_{i=1}^n1 \sum_{i'=1}^n1 h(X_i, X_i') 
  #                         + \sum_{i=1}^n1 \sum_{j=1}^n1 h(X_i, Y_j)
  # = \sum_{i=1}^n1 \sum_{j=1}^n1 h(X_i, Y_j)   
  #                              (the first term is zero since h(x,y) = -h(y,x))
  # = TS (say)
  # Similar result holds for group 2. Therefore, 
  # KW = n1 * \bar{R}^(1)^t Sigma_hat^{-1} \bar{R}^(1) + 
  #             n2 * \bar{R}^(2)^t Sigma_hat^{-1} \bar{R}^(2)
  #    = (n1 / (n1^2 * n^2) + n2 / (n2^2 * n^2) ) TS^t Sigma_hat^{-1} TS
  #    = (1 / (n * n1 * n2)) TS^t Sigma_hat^{-1} TS
  
  TS = colSums(colSums(h(dataX, dataY)))
  KW = sum(TS * solve(Sigma_hat, TS)) / (n * n1 * n2)
  
  return(1-pchisq(KW, p))
}


crossMMD2sample <- function(K, n, m) {
  # Computes the p-value for the studentized cross-MMD test from the paper:
  #          https://arxiv.org/pdf/2211.14908.pdf 
  #          Python code implementation for the paper is available at:
  #          https://github.com/sshekhar17/PermFreeMMD
  # Input: 
  #   K: (n+m) * (n+m) matrix whose i,j th entry is k(Z_i, Z_j) where 
  #      k is the kernel and Z_i are the pooled X and Y observations
  #   n: number of X observations
  #   m: number of Y observations
  # Returns: 
  #   p-value of the test
  
  
  # split the dataset into two equal parts
  n1 = n %/% 2
  m1 = m %/% 2
  
  # compute the Gram matrices
  Kxx = K[1:n1, (n1+1):n, drop = FALSE]
  Kyy = K[(n+1):(n+m1), (n+m1+1):(n+m), drop = FALSE]
  Kxy = K[1:n1, (n+m1+1):(n+m), drop = FALSE]
  Kyx = K[(n+1):(n+m1), (n1+1):n, drop = FALSE]
  
  # compute the numerator of the statistic (main component of the test 
  # statistic)
  Ux = mean(Kxx) - mean(Kxy)
  Uy = mean(Kyx) - mean(Kyy)
  
  U = Ux - Uy
  
  # compute the denominator (required for studentization)
  term1 = (rowMeans(Kxx) - rowMeans(Kxy) - Ux)^2
  sigX2 = mean(term1)
  
  term2 = (rowMeans(Kyx) - rowMeans(Kyy) - Uy)^2
  sigY2 = mean(term2)
  
  sig = sqrt(sigX2 / n1 + sigY2 / m1)
  
  if (!(sig > 0)) {
    print(list(
      term1 = term1,
      term2 = term2,
      sigX2 = sigX2,
      sigY2 = sigY2
    ))
    
    stop(paste("The denominator is", sig))
  }
  
  # obtain the statistic
  T <- U / sig
  return(1-pnorm(T))
}



median_bandwidth <- function(D) {
  # Median bandwidth selection for kernel mmd based tests
  #
  # Input:
  #    D: distance matrix of pooled sample
  # Returns: 
  #    median bandwidth sigma for kernel mmd based tests
  
  dvals = D[upper.tri(D)]
  return(median(dvals))
}




