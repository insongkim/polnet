
#Choosing to not include N_row and N_col as parameters to R function

#'@param edges Matrix of connection strengths as counts
#'@param D The dimensionality of the latent space, 2 dimensions is recommended
#'@return A trained stanmodel object

#'@import Rcpp
#'@import methods
#'@import rstantools
#'@import rstan
#'@useDynLib polnet, .registration = TRUE
#'@export
LSNM <- function(edges, D = 2, ...){
  #model frame
  edge_mat = as.matrix(edges)
  #parameters necessary to run lmem.stan function
  stanlist <- list(edges = edge_mat, D = D, N_row = nrow(edge_mat), N_col = ncol(edge_mat))
  sample_post <- rstan::sampling(stanmodels$LSNM, data = stanlist, ...)
  out <- list(stan_fitted_model = sample_post)
  class(out) <- 'LSNM'
  return(out)
}
