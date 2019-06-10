#' @name LSNM
#' @param edges Matrix or data.frame or igraph of connection strengths as counts
#' @param D The dimensionality of the latent space, 2 dimensions is recommended
#' @param method One of \code{vi} (variational inference) or
#' \code{mcmc} specifying the method of inference. The default is
#' \code{vi}.
#' @param group1.id A character string indicating the name of group1 identifier
#' variable in the \code{edges} data.frame. It is required in the case of data.frame.
#' @param group2.id A character string indicating the name of group2 identifier
#' variable in the \code{edges} data.frame. It is required in the case of data.frame.
#' @param count.id A character string indicating the name of count identifier
#' variable in the \code{edges} data.frame. The variable must be numeric.
#' @param N_fixed_row
#' @param N_fixed_col
#' @param fixed_row_index
#' @param fixed_row_embedding
#' @param fixed_col_index
#' @param fixed_col_embedding
#' @return A trained stanmodel object

#' @import Rcpp
#' @import methods
#' @import rstantools
#' @import rstan
#' @import igraph
#' @import tidyr
#' @useDynLib polnet, .registration = TRUE
#' @export LSNM
#' @examples \dontrun{
#' set.seed(11)
#' sim.data <- random_LSNM_data_cluster(n.cluster=4, group1.center=rbind(c(-0.5,-1), c(-1, 0.3), c(0.4, 1), c(0.2, -0.2))*5, group2.center=rbind(c(-0.5,-1), c(-1, 0.3), c(0.4, 1), c(0.2, -0.2))*5, v=3, sigma_sq_L = 0.5, sigma_sq_P = 0.7, tau=c(0.5, 0.8))
#' res <- LSNM(sim.data$LSNM_data$A, D=2, method = "vi", iter=50000)
#' plot.compare.LSNM(res, sim.data$LSNM_data$Theta, sim.data$LSNM_data$Psi, sim.data$group1.popularity, sim.data$group2.popularity, sim.data$group1.cluster, sim.data$group2.cluster, legend_position = "center")
#' }
#'
LSNM <- function(edges,
                 D = 2,
                 link_function = "poisson",
                 n = NULL,
                 method = c("vi", "mcmc"),
                 group1.id = NULL,
                 group2.id = NULL,
                 count.id = NULL,
                 N_fixed_row=0,
                 N_fixed_col=0,
                 fixed_row_index=vector(),
                 fixed_row_embedding=matrix(0, nrow=0, ncol=2),
                 fixed_col_index=vector(),
                 fixed_col_embedding=matrix(0, nrow=2, ncol=0),
                 ...){

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
  if (!link_function %in% c("poisson", "binomial", "bernoulli"))
    stop("Invalid link function")
  if (link_function=="binomial"&is.null(n)) {
    stop("'n' should be provided")
  }
  if (!method %in% c("vi","mcmc"))
    stop("'method' should be either 'vi' or 'mcmc'")

  if (link_function=="bernoulli") n <- 1
  
  # Input data
  if (class(edges)=="matrix") {
    edge_mat <- edges
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
  } 
  
  if (link_function == "poisson") {
    if (method == "vi") {
      # Parameters necessary to run stan function
      stanlist <- list(edges = edge_mat, D = D, N_row = nrow(edge_mat), N_col = ncol(edge_mat),
                       N_fixed_row=N_fixed_row, N_fixed_col=N_fixed_col, fixed_row_index=fixed_row_index,
                       fixed_row_embedding=fixed_row_embedding, fixed_col_index=fixed_col_index,
                       fixed_col_embedding=fixed_col_embedding)
      sample_post <- rstan::vb(stanmodels$LSNM, data = stanlist, ...)
    } else {
      # Parameters necessary to run stan function
      stanlist <- list(edges = edge_mat, D = D, N_row = nrow(edge_mat), N_col = ncol(edge_mat),
                       N_fixed_row=N_fixed_row, N_fixed_col=N_fixed_col, fixed_row_index=fixed_row_index,
                       fixed_row_embedding=fixed_row_embedding, fixed_col_index=fixed_col_index,
                       fixed_col_embedding=fixed_col_embedding)
      sample_post <- rstan::sampling(stanmodels$LSNM, data = stanlist, ...)
    }
  } else {
    if (method == "vi") {
      # Parameters necessary to run stan function
      stanlist <- list(edges = edge_mat, D = D, n = n, N_row = nrow(edge_mat), N_col = ncol(edge_mat),
                       N_fixed_row=N_fixed_row, N_fixed_col=N_fixed_col, fixed_row_index=fixed_row_index,
                       fixed_row_embedding=fixed_row_embedding, fixed_col_index=fixed_col_index,
                       fixed_col_embedding=fixed_col_embedding)
      sample_post <- rstan::vb(stanmodels$LSNMbinom, data = stanlist, ...)
    } else {
      # Parameters necessary to run stan function
      stanlist <- list(edges = edge_mat, D = D, n = n, N_row = nrow(edge_mat), N_col = ncol(edge_mat),
                       N_fixed_row=N_fixed_row, N_fixed_col=N_fixed_col, fixed_row_index=fixed_row_index,
                       fixed_row_embedding=fixed_row_embedding, fixed_col_index=fixed_col_index,
                       fixed_col_embedding=fixed_col_embedding)
      sample_post <- rstan::sampling(stanmodels$LSNMbinom, data = stanlist, ...)
    }
  }
  

  out <- list(stan_fitted_model = sample_post)

  class(out) <- 'LSNM'
  return(out)
}

#' Get summaries of a LSNM object
#'
#' \code{summary.LSNM()} takes an object returned by
#' \code{LSNM}, and returns a matrix of the mean,
#' standard deviation, and credible interval of the latent space with all chains being merged
#'
#' @param LSNM_Object A trained object of class LSNM
#' @param low_perc The bottom range of the desired credible interval, defaults to 0.1
#' @param high_perc The top range of the credible interval, defaults to 0.9
#' @return A matrix that includes the mean, standard deviation, and credible interval of the latent space estimated by the LSNM algorithm. The row embeddings are the client latent space positions, while the column embeddings are the legislator latent space positions.
#'

#' @useDynLib polnet, .registration = TRUE
#' @export

summary.LSNM <- function(LSNM_Object,
                         low_perc = 0.1,
                         high_perc = 0.9){

  df_fit <- as.data.frame(LSNM_Object$stan_fitted_model)
  nms <- df_fit[ , grepl( "^col_embedding|^row_embedding" , names(df_fit) )]
  final_df <- as.data.frame(cbind(colMeans(nms), apply(nms, 2, sd), apply(nms, 2, quantile, low_perc), apply(nms, 2, quantile, high_perc)))
  names(final_df) <- c('Mean', 'SD', '10%', '90%')

  return(final_df)
}


#' Plot the posterior means of a LSNM object
#'
#' \code{plot.LSNM()} takes an object returned by
#' \code{LSNM}, and returns a plot of the posterior means.
#' Use standard arguments to the \code{plot} function to modify the plot as needed.
#'
#' @param LSNM_Object A trained object of class LSNM
#' @param group1_cluster A vector representing the cluster of group1
#' @param group2_cluster A vector representing the cluster of group2
#' @return a plot of the posterior means
#' @useDynLib polnet, .registration = TRUE
#' @export

plot.LSNM <- function(LSNM_Object,
                      group1_cluster = NULL,
                      group2_cluster = NULL,
                      main = "Estimated LSNM Positions",
                      legend = c("Group1", "Group2"),
                      legend_position = "topleft",
                      ...){

  m <- LSNM_Object$stan_fitted_model@par_dims$row_factor_adj # number of group1
  n <- LSNM_Object$stan_fitted_model@par_dims$col_factor_adj # number of group2
  D <- LSNM_Object$stan_fitted_model@par_dims$cov_embedding_diag # number of dimensions

  if (is.null(group1_cluster)) group1_cluster <- rep("black", m)
  if (is.null(group2_cluster)) group2_cluster <- rep("black", n)

  df_fit <- as.data.frame(LSNM_Object$stan_fitted_model)
  nms <- df_fit[ , grepl( "^col_embedding|^row_embedding|^col_factor|^row_factor" , names(df_fit) )]
  plot.data <- colMeans(nms) # posterior mean

  if (D==1) {
    row_size <- exp(plot.data[paste0("row_factor_adj[",1:m,"]")]) # size of group1
    row_size <- 2*row_size/max(row_size)
    col_size <- exp(plot.data[paste0("col_factor_adj[",1:n,"]")]) # size of group2
    col_size <- 2*col_size/max(col_size)

    row_elements <- plot.data[paste0("row_embedding[",1:m,",1]")]
    col_elements <- plot.data[paste0("col_embedding[1,",1:n,"]")]


    positions <- sort(c(row_elements,col_elements))
    row_ord <- which(substr(names(positions),1,3)=="row")
    col_ord <- which(substr(names(positions),1,3)=="col")

    row_pos <- positions[row_ord]
    col_pos <- positions[col_ord]

    row_lab <- paste0(substr(legend[1],1,1), order(row_elements))
    col_lab <- paste0(substr(legend[2],1,1), order(col_elements))

    plot(x = row_pos,
         y = row_ord,
         pch = 1,
         cex = row_size,
         xlab = "Latent Space Dimension 1",
         ylab = "",
         yaxt = "n",
         main = main,
         col = group1_cluster,...)
    axis(side = 2,
         at = row_ord,
         labels = row_lab)
    points(x = col_pos,
           y = col_ord,
           pch = 0,
           cex = col_size,
           col = group2_cluster)
    axis(side = 4,
         at = col_ord,
         labels = col_lab)
    legend(legend_position,
           legend = legend,
           pch = c(1, 0),
           cex = 0.8,
           box.lty = 0,
           inset = 0.01)


  } else {
    row_size <- exp(plot.data[paste0("row_factor_adj[",1:m,"]")]) # size of group1
    row_size <- 2*row_size/max(row_size)
    col_size <- exp(plot.data[paste0("col_factor_adj[",1:n,"]")]) # size of group2
    col_size <- 2*col_size/max(col_size)

    plot(x = plot.data[paste0("row_embedding[",1:m,",1]")],
         y = plot.data[paste0("row_embedding[",1:m,",2]")],
         pch = 1,
         cex = row_size,
         xlab = "Latent Space Dimension 1",
         ylab = "Latent Space Dimension 2",
         main = main,
         col = group1_cluster,...)
    points(x = plot.data[paste0("col_embedding[1,",1:n,"]")],
           y = plot.data[paste0("col_embedding[2,",1:n,"]")],
           pch = 0,
           cex = col_size,
           col = group2_cluster)
    legend(legend_position,
           legend = legend,
           bg = "transparent",
           pch = c(1, 0),
           cex = 0.8,
           box.lty = 0,
           inset = 0.01)
  }

}
