#' @name LSNM
#' @param edges Matrix or data.frame of connection strengths as counts
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
#' @importFrom rstan vb
#' @importFrom rstan sampling
#' @useDynLib polnet, .registration = TRUE
#' @export
LSNM <- function(edges,
                 D = 2,
                 method = c("vi", "mcmc"),
                 group1.id = NULL,
                 group2.id = NULL,
                 count.id = NULL,
                 N_fixed_row=0,
                 N_fixed_col=0,
                 fixed_row_index=NULL,
                 fixed_row_embedding=NULL,
                 fixed_col_index=NULL,
                 fixed_col_embedding=NULL,
                 ...){

  ## Warning for missing parameter
  if (missing(edges))
    stop("'edges' should be provided")
  if (!class(edges)%in%c("matrix","data.frame"))
    stop("'edges' should be either matrix or data.frame")
  if (class(edges)=="data.frame"&is.null(group1.id))
    stop("'group1.id' should be provided")
  if (class(edges)=="data.frame"&is.null(group2.id))
    stop("'group2.id' should be provided")
  if (class(edges)=="data.frame"&is.null(count.id))
    stop("'count.id' should be provided")
  if (!method %in% c("vi","mcmc"))
    stop("'method' should be either 'vi' or 'mcmc'")

  # Input data
  if (class(edges)=="matrix") {
    edge_mat <- edges
  } else {
    edge_mat <- matrix(data = edge$count.id, nrow = nrow(edges), ncol = ncol(edges))
    rownames(edge_mat) <- edges$group1
    colnames(edge_mat) <- edges$group2
  }

  if (method == "vi") {
    # Parameters necessary to run stan function
    stanlist <- list(edges = edge_mat, D = D, N_row = nrow(edge_mat), N_col = ncol(edge_mat))
    sample_post <- rstan::vb(stanmodels$LSNMshort, data = stanlist, ...)
  } else {
    # Parameters necessary to run stan function
    stanlist <- list(edges = edge_mat, D = D, N_row = nrow(edge_mat), N_col = ncol(edge_mat),
                     N_fixed_row=N_fixed_row, N_fixed_col=N_fixed_col, fixed_row_index=fixed_row_index,
                     fixed_row_embedding=fixed_row_embedding, fixed_col_index=fixed_col_index,
                     fixed_col_embedding=fixed_col_embedding)
    sample_post <- rstan::sampling(stanmodels$LSNM, data = stanlist, ...)
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
#' @return a plot of the posterior means
#' @useDynLib polnet, .registration = TRUE
#' @export

plot.LSNM <- function(LSNM_Object,
                      main = "Estimated LSNM Positions",
                      legend = c("Group1", "Group2"),
                      legend_position = "topleft",
                      ...){

  m <- LSNM_Object$stan_fitted_model@par_dims$row_factor # number of group1
  n <- LSNM_Object$stan_fitted_model@par_dims$col_factor # number of group2
  D <- LSNM_Object$stan_fitted_model@par_dims$mu_col_embedding # number of dimensions

  df_fit <- as.data.frame(LSNM_Object$stan_fitted_model)
  nms <- df_fit[ , grepl( "^col_embedding|^row_embedding|^col_factor|^row_factor" , names(df_fit) )]
  plot.data <- colMeans(nms) # posterior mean

  if (D==1) {
    row_size <- exp(plot.data[paste0("row_factor[",1:m,"]")]) # size of group1
    row_size <- 2*row_size/max(row_size)
    col_size <- exp(plot.data[paste0("col_factor[",1:n,"]")]) # size of group2
    col_size <- 2*col_size/max(col_size)

    row_elements <- plot.data[paste0("row_embedding[",1:m,",1]")]
    col_elements <- plot.data[paste0("col_embedding[",1:n,",1]")]


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
         main = main, ...)
    axis(side = 2,
         at = row_ord,
         labels = row_lab)
    points(x = col_pos,
           y = col_ord,
           pch = 16,
           cex = col_size,
           col = rgb(0,0,0,alpha=0.8))
    axis(side = 4,
         at = col_ord,
         labels = col_lab)
    legend(legend_position,
           legend = legend,
           pch = c(1, 16),
           col = c("black", rgb(0,0,0,alpha=0.8)),
           cex = 0.8,
           box.lty = 0,
           inset = 0.01)


  } else {
    row_size <- exp(plot.data[paste0("row_factor[",1:m,"]")]) # size of group1
    row_size <- 2*row_size/max(row_size)
    col_size <- exp(plot.data[paste0("col_factor[",1:n,"]")]) # size of group2
    col_size <- 2*col_size/max(col_size)

    plot(x = plot.data[paste0("row_embedding[",1:m,",1]")],
         y = plot.data[paste0("row_embedding[",1:m,",2]")],
         pch = 1,
         cex = row_size,
         xlab = "Latent Space Dimension 1",
         ylab = "Latent Space Dimension 2",
         main = main, ...)
    points(x = plot.data[paste0("col_embedding[",1:n,",1]")],
           y = plot.data[paste0("col_embedding[",1:n,",2]")],
           pch = 16,
           cex = col_size,
           col = rgb(0,0,0,alpha=0.8))
    legend(legend_position,
           legend = legend,
           pch = c(1, 16),
           col = c("black", rgb(0,0,0,alpha=0.8)),
           cex = 0.8,
           box.lty = 0,
           inset = 0.01)
  }

}
