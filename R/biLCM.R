#'@name biLCM
#'@param A Matrix of connection strength as counts
#'@param m Number of clients
#'@param n Number of politicians
#'@param k Number of link communities
#'@return A list specifying the membership weights, client membership distributions, and politician membership distributions

#'@export

biLCM <- function(edges, 
                  group1.id = NULL,
                  group2.id = NULL,
                  count.id = NULL
                  k = NULL,
                  seed = 12345,
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
  set.seed(seed)
  kappa <- rgamma(k, 2, 0.5)
  alpha <- matrix(runif(m*k), m, k)
  alpha[,] <- alpha/rowSums(alpha)
  beta <- matrix(runif(n*k), n, k)
  beta[,] <- beta/rowSums(beta)
  q <- array(0, c(m, n, k))
  
  iter_count <- 0
  llik_iter <- numeric(max.iter)
  
  repeat {
    iter_count <- iter_count + 1
    # E-step

    
    # M-step

    
    # Calculate log likelihood
    llik_new <- 0
    llik_iter[iter_count] <- llik_new
    
    ## check convergence
    if (llik_new - llik_old < tolerance || iter_count == max.iter) break
    ## update log-likelihood
    llik_old <- llik_new
  }
    
	answer <- list(kappa = kappa, alpha = alpha, beta = beta)
	return(answer)
}
