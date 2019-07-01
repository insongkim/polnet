#'@param m The number of group1
#'@param n The number of group2
#'@param k The number of link communities
#'@param kappa_weight The base weight. Expecting a k length vector. Assumed to be drawn from a uniform distribution with lower limits of kappa_lb and upper limits of keppa_ub. An optional argument, but must specifiy this or kappa_lb and kappa_ub.
#'@param kappa_lb The lower limits of kappa. An optional argument, required if missing kappa.
#'@param kappa_ub The upper limits of kappa. An optional argument, required if missing kappa.
#'@param alpha_membership The desired latent mixed-membership of the group1. Expecting a matrix with m rows and k columns. An optional argument, but must specify this or alpha_num.
#'@param beta_membership The desired latent mixed-membership of the group2. Expecting a matrix with n rows and k columns. An optional argument, but must specify this or beta_num.
#'@param alpha_num The desired number of communities that each member of group1 is affiliated with. Assumed to be drawn from a poisson distribution with mean of k/2. Note that 0 is replaced with 1 and the value that exceeds k is replaced with k. An optional argument, required if missing alpha_membership.
#'@param beta_num The desired number of communities that each member of group2 is affiliated with. Assumed to be drawn from a poisson distribution with mean of k/2. Note that 0 is replaced with 1 and the value that exceeds k is replaced with k. An optional argument, required if missing beta_membership.
#'@return A list of kappa_weight, alpha_membership, beta_membership, and then the randomly generated poisson count matrix A.

#'@export

random_biLCM_data <- function(m, 
                              n, 
                              k,
                              kappa_weight = NULL,
                              kappa_lb = NULL,
                              kappa_ub = NULL,
                              alpha_membership = NULL,
                              beta_membership = NULL,
                              alpha_num = NULL,
                              beta_num = NULL) {
  
  # Generating kappa
  if(is.null(kappa_weight)){
    if(is.null(kappa_lb) | is.null(kappa_ub)){
      stop("Invalid Kappa Parameters")
    }else{
      kappa_weight <- runif(k, min = kappa_lb, max = kappa_ub)
    }
  }

  # Generating alpha
  if(is.null(alpha_membership)){
    if(is.null(alpha_num)){
      alpha_num <- rpois(m, k/2)
      alpha_num[which(alpha_num==0)] <- 1
      alpha_num[which(alpha_num>k)] <- k
    }
    alpha_membership <- matrix(0, m, k)
    alpha_index <- sapply(alpha_num, function(x) sample(1:k, x, replace = F))
    for (i in 1:length(alpha_index)) {
      alpha_membership[i, alpha_index[[i]]] <- runif(length(alpha_index[[i]]), min = 0.01)
    }    
    alpha_membership <- alpha_membership/rowSums(alpha_membership)
  }
  
  # Generating beta
  if(is.null(beta_membership)){
    if(is.null(beta_num)){
      beta_num <- rpois(n, k/2)
      beta_num[which(beta_num==0)] <- 1
      beta_num[which(beta_num>k)] <- k
    }
    beta_membership <- matrix(0, n, k)
    beta_index <- sapply(beta_num, function(x) sample(1:k, x, replace = F))
    for (i in 1:length(beta_index)) {
      beta_membership[i, beta_index[[i]]] <- runif(length(beta_index[[i]]), min = 0.01)
    }
    beta_membership <- beta_membership/rowSums(beta_membership)
  }
  
  mu <- sweep(alpha_membership, 2, kappa_weight, "*") %*% t(beta_membership)
  A_mat <- matrix(sapply(mu, function(x) rpois(1, x)), m, n)
  
  return(list(Kappa = kappa_weight, Alpha = alpha_membership, Beta = beta_membership, A = A_mat))
}
