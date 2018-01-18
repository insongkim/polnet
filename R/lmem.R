#'Estimates regression params using a normal model for errors, 
#'a normal model for Beta, and a Cauchy for sigma

#'@param formula The formula for the regression. Dependent variable ~ independent variables
#'@param data The data used for the regression, data frame
#'@param b_loc Prior for Beta mean
#'@param b_scale Prior for Beta sd
#'@param sigma_scale Prior for sigma sd (mean is 0)
#'@return A trained stanmodel object

#'@import Rcpp
#'@import methods
#'@import rstantools
#'@import rstan
#'@useDynLib lmem, .registration = TRUE
#'@export
lmem <- function(formula, data,   b_loc = 0, b_scale = 1000, sigma_scale = 50, ...){
  #model frame
  frame <- lm(formula, data, method = 'model.frame')
  #construct design matrix
  X <- model.matrix(formula, data = frame)
  #separate X and Y
  if (colnames(X)[1] == "(Intercept") X <- X[, -1, drop = FALSE]
  Y <- model.response(frame, type = "numeric")
  #parameters necessary to run lmem.stan function
  stanlist <- list(n = nrow(X), y = Y, k = ncol(X), X = X, b_loc = b_loc, b_scale = b_scale, sigma_scale = sigma_scale)
  sample_post <- rstan::sampling(stanmodels$lmem, data = stanlist, ...)
  out <- list(stan_fitted_model = sample_post)
  class(out) <- 'lmem'
  return(out)
}