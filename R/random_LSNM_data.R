#'@param D The dimensionality of the latent space, 2 dimensions is recommended
#'@param clients The number of clients 
#'@param legislators The number of legislators
#'@param client_space The desired latent space positions of the clients. Expecting a matrix with clients rows and D columns. Assumed to be drawn from a multivariate normal distribution with a mean of 0. An optional argument, but must specify this or tau. 
#'@param legislator_space The desired latent space positions of the legislators. Expecting a matrix with legislators rows and D columns. Assumed to be drawn from a multivariate normal distribution. An optional argument, but must specifiy this or mu and Sigma. 
#'@param mu The mean vector to generate legislator_space. Expecting a vector of length D. An optional argument, required if missing legislator_space.
#'@param Sigma The variance-covariance matrix used to generate legislator_space. Expecting a D by D matrix with a reasonable variance-covariance structure. An optional argument, required if missing legislator_space.
#'@param tau The variance-covariance matrix used to generate client_space. Expecting a D by D matrix with a diagonal variance-covariance structure. An optional argument, required if missing client_space. 
#'@param alpha_popularity The client popularity factors used to account for baseline likelihood to lobby. Expecting a clients length vector. Assumed to be drawn from a normal distribution. An optional argument, but must specify this or v and sigma_sq_L
#'@param Beta_popularity The legislator popularity factors used to account for baseline likelihood to sponsor bills. Expecting a legislators length vector. Assumed to be drawn from a normal distribution with mean 0. An optional argument, but must specify this or sigma_sq_P
#'@param v The mean used to generate the alpha_popularity vector. An optional argument, required if missing alpha_popularity
#'@param sigma_sq_L The variance used to generate the alpha_popularity vector. An optional argument, required if missing alpha_popularity
#'@param sigma_sq_P The variance used to generate the Beta_popularity vector. An optional argument, required if missing Beta_popularity
#'@return A list of the latent client space, the latent legislator space, and then the randomly generated poisson count matrix A. 

#'@import MASS
#'@useDynLib polnet, .registration = TRUE
#'@export


random_LSNM_data <- function(D, clients, legislators, client_space = NULL, legislator_space = NULL, mu = NULL, Sigma = NULL, tau = NULL, alpha_popularity = NULL, Beta_popularity = NULL, v = NULL, sigma_sq_L = NULL, sigma_sq_P = NULL){
  #Generate the alpha, Beta vectors if not existing already
  if(is.null(alpha_popularity)){
    if(is.null(v) | is.null(sigma_sq_L)){
      stop("Invalid Alpha Parameters")
    }else{
      alpha_popularity = rnorm(clients, v, sigma_sq_L)
    }
  }
  if(is.null(Beta_popularity)){
    if(is.null(sigma_sq_P)){
      stop("Invalid Beta Parameters")
    }else{
      Beta_popularity = rnorm(legislators, 0, sigma_sq_P)
    }
  }
  #putting together matrix of alpha + Beta
  a_B = outer(alpha_popularity, Beta_popularity, '+')
  if(is.null(client_space)){
    if(is.null(tau)){
      stop('Invalid Theta Parameters')
    }else{
      if(!all(tau[lower.tri(tau)] == 0, tau[upper.tri(tau)] == 0)){
        stop('Tau covariance matrixmust be diagonal')
      }
      client_space = mvrnorm(n = clients, mu = rep(0, D), Sigma = tau)
    }
  }
  if(is.null(legislator_space)){
    if(is.null(mu) | is.null(Sigma)){
      stop('Invalid Psi Parameters')
    }else{
      legislator_space = mvrnorm(n = legislators, mu = mu, Sigma = Sigma)
    }
  }
  distances = outer(split(client_space, row(client_space)), split(legislator_space, row(legislator_space)), Vectorize(l1_norm))
  diff_vec = a_B - distances
  A_mat = matrix(sapply(as.vector(diff_vec), function(x) rpois(1, exp(x))), nrow = nrow(distances))
  return(list(Theta = client_space, Psi = legislator_space, A = A_mat))
}


l1_norm <- function(vec1, vec2){
  return(sum(abs(vec1 - vec2)))
}
