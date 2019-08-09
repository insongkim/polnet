#'@name biLCM
#'@param edges Matrix or data.frame or igraph of connection strengths as counts (NA is considered as no edges)
#'@param group1.id A character string indicating the name of group1 identifier
#'variable in the \code{edges} data.frame. It is required in the case of data.frame.
#'@param group2.id A character string indicating the name of group2 identifier
#'variable in the \code{edges} data.frame. It is required in the case of data.frame.
#'@param count.id A character string indicating the name of count identifier
#'variable in the \code{edges} data.frame. The variable must be numeric.
#'@param k Number of link communities.
#'@param tolerance Tolerance for convergence.
#'@param max.iter Maximum value of iteration numbers. It is required to avoid infinite loop.
#'@return A list specifying the log-likelihood for each iteration, membership weights, membership distributions of group1 and group2

#'@export
#'@examples \dontrun{
#'set.seed(11)
#'sim.data <- random_biLCM_data(m = 100, n = 50, k = 4, kappa_weight = NULL, a = 10000, b = 1, alpha_membership = NULL, alpha_c = rep(0.5, 100), beta_membership = NULL, beta_c = rep(0.5, 50), non_zero = TRUE)
#'res <- biLCM(edges = sim.data$A, group1.id = NULL, group2.id = NULL, count.id = NULL, k = 4, tolerance = 1e-6, max.iter = 200)
#'plot.compare.biLCM(res, sim.data, group1 = TRUE, nth = 1)
#' }

biLCM <- function(edges, 
                  group1.id = NULL,
                  group2.id = NULL,
                  count.id = NULL,
                  k = NULL,
                  tolerance = 1e-6,
                  max.iter = 200){
  
  ## Warning for missing parameter
  if (missing(edges))
    stop("'edges' should be provided")
  if (!class(edges)%in%c("matrix","data.frame", "igraph"))
    stop("'edges' should be matrix or data.frame or igraph")
  if (class(edges)=="data.frame"&is.null(group1.id))
    stop("'group1.id' should be provided")
  if (class(edges)=="data.frame"&is.null(group2.id))
    stop("'group2.id' should be provided")
  if (class(edges)%in%c("data.frame", "igraph")&is.null(count.id))
    stop("'count.id' should be provided")
  if (is.null(k))
    stop("'k' should be provided")
  
  # Input data
  if (class(edges)=="matrix") {
    edge_mat <- edges
    edge_mat[is.na(edge_mat)] <- 0
  } else {
    
    if (class(edges)=="igraph") {
      edges <- igraph::as_data_frame(edges, what = "edges")
      group1.id <- "from"
      group2.id <- "to"
    }
    
    edges <- edges[,c(group1.id, group2.id, count.id)]
    edge_mat <- tidyr::spread(edges, group2.id, count.id)
    rownames(edge_mat) <- edge_mat[,group1.id]
    edge_mat <- edge_mat[,-1]
    edge_mat <- as.matrix(edge_mat)
    edge_mat[is.na(edge_mat)] <- 0
  } 
  
  m <- nrow(edge_mat)
  n <- ncol(edge_mat)
  
  # Initialize parameters
  kappa <- rgamma(k, 2, 0.5)
  alpha <- t(rdirichlet(k, rep(1, m)))
  beta <- t(rdirichlet(k, rep(1, n)))
  q <- array(0, c(m, n, k))
  
  iter_count <- 0
  llik_iter <- numeric(max.iter)
  
  kab <- sweep(alpha, 2, kappa, "*") %*% t(beta)
  llik_old <- sum(edge_mat * log(kab)) - sum(kab)
  
  repeat {
    iter_count <- iter_count + 1
    
    ## E-step
    numer <- sapply(1:k, function(z) kappa[z]*outer(alpha[,z], beta[,z], "*"), simplify = "array")
    log_numer <- log(numer)
    c <- apply(log_numer, c(1,2), max)
    c <- replicate(k, c)
    numer <- exp(log_numer-c)
    denom <- apply(numer, c(1,2), sum)
    q <- sweep(numer, c(1,2), denom, "/")
    
    ## M-step
    # Update alpha
    numer <- sapply(1:k, function(z) rowSums(sweep(edge_mat, c(1,2), q[,,z], "*")))
    denom <- colSums(numer)
    alpha <- sweep(numer, 2, denom, "/")
    # Update beta
    numer <- sapply(1:k, function(z) colSums(sweep(edge_mat, c(1,2), q[,,z], "*")))
    beta <- sweep(numer, 2, denom, "/")
    # Update kappa
    numer <- denom
    kappa <- numer/sapply(1:k, function(z) sum(outer(alpha[,z], beta[,z], "*")))
    
    # Calculate log likelihood
    kab <- sweep(alpha, 2, kappa, "*") %*% t(beta)
    log_kab <- log(kab)
    log_kab[which(is.na(log_kab))] <- 0
    log_kab[which(is.infinite(log_kab))] <- -1.797693e+308
    llik_new <- sum(edge_mat * log_kab) - sum(kab)
    llik_iter[iter_count] <- llik_new
    
    # Check convergence
    if (llik_new - llik_old < tolerance || iter_count == max.iter) break
    # Update log-likelihood
    llik_old <- llik_new
  }
  
  out <- list(loglikelihood = llik_iter, kappa = kappa, alpha = alpha, beta = beta)
  
  class(out) <- "biLCM"
  return(out)
}
