#' @name LSNM
#' @param edges Matrix or data.frame or igraph of connection strengths as counts (NA is considered as no edges)
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
                 method = c("vi", "mcmc", "vi-mcmc"),
                 group1.id = NULL,
                 group2.id = NULL,
                 count.id = NULL,
                 N_fixed_row=0,
                 N_fixed_col=0,
                 fixed_row_index=vector(),
                 fixed_row_embedding=matrix(0, nrow=0, ncol=2),
                 fixed_col_index=vector(),
                 fixed_col_embedding=matrix(0, nrow=2, ncol=0),
                 fixed.actor.object=NULL,
                 iter.vb=NULL, iter.mcmc=NULL,
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
  if (!method %in% c("vi","mcmc","vi-mcmc"))
    stop("'method' should be either 'vi' or 'mcmc' or 'vi-mcmc")
  if (!is.null(fixed.actor.object) & class(fixed.actor.object)!="LSNM_fixed_actors")
    stop("fixed.actor.object must be of class LSNM_fixed_actors.")
  if (missing(iter) & is.null(iter.vb) & is.null(iter.mcmc)){
    iter <- 2000
    warning("Number of iteration=2000 applies by default.")
  }
  if (is.null(iter.vb))
    iter.vb <- iter
  if (is.null(iter.mc.c))
    iter.mcmc <- iter

  if (link_function=="bernoulli") n <- 1
  
  ## Create fixed.actor.list if NULL
  if (is.null(fixed.actor.object)){
    fixed.actor.object <- list(N_fixed_row=N_fixed_row,
                               N_fixed_col=N_fixed_col,
                               fixed_row_index=fixed_row_index,
                               fixed_row_embedding=fixed_row_embedding,
                               fixed_col_index=fixed_col_index,
                               fixed_col_embedding=fixed_col_embedding)
    class(fixed.actor.object) <- "LSNM_fixed_actors"
  }
  

  
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
  
  if (link_function == "poisson") {
    if (method == "vi") {
      # Parameters necessary to run stan function
      stanlist <- list(edges = edge_mat, D = D, N_row = nrow(edge_mat), N_col = ncol(edge_mat),
                       N_fixed_row=fixed.actor.object$N_fixed_row, N_fixed_col=fixed.actor.object$N_fixed_col, fixed_row_index=fixed.actor.object$fixed_row_index,
                       fixed_row_embedding=fixed.actor.object$fixed_row_embedding, fixed_col_index=fixed.actor.object$fixed_col_index,
                       fixed_col_embedding=fixed.actor.object$fixed_col_embedding)
      sample_post <- rstan::vb(stanmodels$LSNM, data = stanlist, iter=iter.vb, ...)
    } else if (method == "mcmc"){
      # Parameters necessary to run stan function
      stanlist <- list(edges = edge_mat, D = D, N_row = nrow(edge_mat), N_col = ncol(edge_mat),
                       N_fixed_row=fixed.actor.object$N_fixed_row, N_fixed_col=fixed.actor.object$N_fixed_col, fixed_row_index=fixed.actor.object$fixed_row_index,
                       fixed_row_embedding=fixed.actor.object$fixed_row_embedding, fixed_col_index=fixed.actor.object$fixed_col_index,
                       fixed_col_embedding=fixed.actor.object$fixed_col_embedding)
      sample_post <- rstan::sampling(stanmodels$LSNM, data = stanlist, iter=iter.mcmc, ...)
    } else if (method == "vi-mcmc"){
      # Parameters necessary to run stan function
      stanlist <- list(edges = edge_mat, D = D, N_row = nrow(edge_mat), N_col = ncol(edge_mat),
                       N_fixed_row=fixed.actor.object$N_fixed_row, N_fixed_col=fixed.actor.object$N_fixed_col, fixed_row_index=fixed.actor.object$fixed_row_index,
                       fixed_row_embedding=fixed.actor.object$fixed_row_embedding, fixed_col_index=fixed.actor.object$fixed_col_index,
                       fixed_col_embedding=fixed.actor.object$fixed_col_embedding)
      sample_post <- rstan::vb(stanmodels$LSNM, data = stanlist, iter=iter.vb, ...)
      
      mcmc.input <- list(stan_fitted_model = sample_post)
      class(mcmc.input) <- 'LSNM'
      
      fixed.actor.object <- choose.fix(mcmc.input, n.wild=1, choose.method="axis")
      
      stanlist <- list(edges = edge_mat, D = D, N_row = nrow(edge_mat), N_col = ncol(edge_mat),
                       N_fixed_row=fixed.actor.object$N_fixed_row, N_fixed_col=fixed.actor.object$N_fixed_col, fixed_row_index=fixed.actor.object$fixed_row_index,
                       fixed_row_embedding=fixed.actor.object$fixed_row_embedding, fixed_col_index=fixed.actor.object$fixed_col_index,
                       fixed_col_embedding=fixed.actor.object$fixed_col_embedding)
      sample_post <- rstan::sampling(stanmodels$LSNM, data = stanlist, iter=iter.mcmc, ...)
    }
  } else {
    if (method == "vi") {
      # Parameters necessary to run stan function
      stanlist <- list(edges = edge_mat, D = D, n = n, N_row = nrow(edge_mat), N_col = ncol(edge_mat),
                       N_fixed_row=fixed.actor.object$N_fixed_row, N_fixed_col=fixed.actor.object$N_fixed_col, fixed_row_index=fixed.actor.object$fixed_row_index,
                       fixed_row_embedding=fixed.actor.object$fixed_row_embedding, fixed_col_index=fixed.actor.object$fixed_col_index,
                       fixed_col_embedding=fixed.actor.object$fixed_col_embedding)
      sample_post <- rstan::vb(stanmodels$LSNMbinom, data = stanlist, iter=iter.vb, ...)
    } else if (method == "mcmc") {
      # Parameters necessary to run stan function
      stanlist <- list(edges = edge_mat, D = D, n = n, N_row = nrow(edge_mat), N_col = ncol(edge_mat),
                       N_fixed_row=fixed.actor.object$N_fixed_row, N_fixed_col=fixed.actor.object$N_fixed_col, fixed_row_index=fixed.actor.object$fixed_row_index,
                       fixed_row_embedding=fixed.actor.object$fixed_row_embedding, fixed_col_index=fixed.actor.object$fixed_col_index,
                       fixed_col_embedding=fixed.actor.object$fixed_col_embedding)
      sample_post <- rstan::sampling(stanmodels$LSNMbinom, data = stanlist, iter=iter.mcmc, ...)
    } else if (method == "vi-mcmc") {
      # Parameters necessary to run stan function
      stanlist <- list(edges = edge_mat, D = D, n = n, N_row = nrow(edge_mat), N_col = ncol(edge_mat),
                       N_fixed_row=fixed.actor.object$N_fixed_row, N_fixed_col=fixed.actor.object$N_fixed_col, fixed_row_index=fixed.actor.object$fixed_row_index,
                       fixed_row_embedding=fixed.actor.object$fixed_row_embedding, fixed_col_index=fixed.actor.object$fixed_col_index,
                       fixed_col_embedding=fixed.actor.object$fixed_col_embedding)
      sample_post <- rstan::vb(stanmodels$LSNMbinom, data = stanlist, iter=iter.vb, ...)
      
      mcmc.input <- list(stan_fitted_model = sample_post)
      class(mcmc.input) <- 'LSNM'
      
      fixed.actor.object <- choose.fix(mcmc.input, n.wild=2, choose.method="axis")
      
      stanlist <- list(edges = edge_mat, D = D, n = n, N_row = nrow(edge_mat), N_col = ncol(edge_mat),
                       N_fixed_row=fixed.actor.object$N_fixed_row, N_fixed_col=fixed.actor.object$N_fixed_col, fixed_row_index=fixed.actor.object$fixed_row_index,
                       fixed_row_embedding=fixed.actor.object$fixed_row_embedding, fixed_col_index=fixed.actor.object$fixed_col_index,
                       fixed_col_embedding=fixed.actor.object$fixed_col_embedding)
      sample_post <- rstan::sampling(stanmodels$LSNMbinom, data = stanlist, iter=iter.mcmc, ...)
      
    }
  }  
  

  out <- list(stan_fitted_model = sample_post)

  class(out) <- 'LSNM'
  return(out)
}
#' Choose the actors whose coordinates will be fixed
#'
#' \code{choose.fixed()} takes an object returned by
#' \code{LSNM}, and returns a list containing the parameters
#' of actors whose coordinates are fixed.
#'
#' @param LSNM_Object A trained object of class LSNM
#' @param n.wild The number of actors to be fixed in each octants
#' @param method The method to choose the wildest actors
#' @return a list containing the parameters of actors whose coordinates are fixed.
#'

#' @useDynLib polnet, .registration = TRUE
#' @export choose.fix
#' 
#' 
choose.fix <- function(LSNM_Object, 
                       n.wild=1, 
                       choose.method="axis"){
    
  df_fit <- as.data.frame(LSNM_Object$stan_fitted_model)
  df_fit.mean <- colMeans(df_fit)
  
  D <- max(as.numeric(gsub("^cov_embedding_diag\\[(\\d+)\\]$", "\\1", colnames(df_fit)[grep("^cov_embedding_diag\\[(\\d+)\\]$", colnames(df_fit))])))
  row.actors <- as.numeric(gsub("^row_embedding\\[(\\d+),1\\]$", "\\1", colnames(df_fit)[grep("^row_embedding\\[\\d+,1\\]$", colnames(df_fit))]))
  col.actors <- as.numeric(gsub("^col_embedding\\[1,(\\d+)\\]$", "\\1", colnames(df_fit)[grep("^col_embedding\\[1,\\d+\\]$", colnames(df_fit))]))
  
  
  if (D==2) {
    row.embeddings <- df_fit.mean[grep("^row_embedding\\[\\d+,\\d+\\]$", names(df_fit.mean))]  
    row.embeddings.df <- cbind(df_fit.mean[paste0("row_embedding[", row.actors, ",1]")],
                               df_fit.mean[paste0("row_embedding[", row.actors, ",2]")])
    rownames(row.embeddings.df) <- 1:nrow(row.embeddings.df)
    row.octants <- apply(row.embeddings.df, 1, find.octants)
    
    col.embeddings <- df_fit.mean[grep("^col_embedding\\[\\d+,\\d+\\]$", names(df_fit.mean))]  
    col.embeddings.df <- cbind(df_fit.mean[paste0("col_embedding[1,", col.actors, "]")],
                               df_fit.mean[paste0("col_embedding[2,", col.actors, "]")])
    rownames(col.embeddings.df) <- 1:nrow(col.embeddings.df)
    col.octants <- apply(col.embeddings.df, 1, find.octants)
    
    if (choose.method=="octant"){
      
      ### Find the Wildest Row Actors, prefarbly in Different Octants
      
      ix1 <- order(apply(row.embeddings.df, 1, min_dist_to_octant_line),decreasing=T)[1]
      octant1 <- row.octants[ix1]
      
      octant1.pool <- octant1
      octant1.found <- FALSE
      
      while (!octant1.found){
        ix1.pool <- as.numeric(names(row.octants[row.octants %in% octant1.pool]))
        ix1 <- as.numeric(names(row.embeddings.df[ix1.pool,1])
                          [order(apply(row.embeddings.df[ix1.pool,], 1, min_dist_to_octant_line), decreasing = TRUE)][1:n.wild])
        ix1 <- ix1[!is.na(ix1)]
        if (length(ix1)>=n.wild){
          octant1 <- row.octants[ix1]
          octant1.found <- TRUE
        } else octant1.pool <- c(octant1.pool, (octant1.pool-1) %% 8 , (octant1.pool+1) %% 8)
      }
      
      octant2.pool <- unique(octant1 + 4)
      octant2.found <- FALSE
      
      while (!octant2.found){
        ix2.pool <- as.numeric(names(row.octants[row.octants %in% octant2.pool]))
        ix2 <- as.numeric(names(row.embeddings.df[ix2.pool,1])
                          [order(apply(row.embeddings.df[ix2.pool,], 1, min_dist_to_octant_line), decreasing = TRUE)][1:n.wild])
        ix2 <- ix2[!is.na(ix2)]
        if (length(ix2)>=n.wild){
          octant2 <- row.octants[ix2]
          octant2.found <- TRUE
        } else octant2.pool <- c(octant2.pool, (octant2.pool-1) %% 8 , (octant2.pool+1) %% 8)
      }
      
      octant3.pool <- unique((octant1+1) %% 8 + 1)
      octant3.found <- FALSE
      
      while (!octant3.found){
        ix3.pool <- as.numeric(names(row.octants[row.octants %in% octant3.pool]))
        ix3 <- as.numeric(names(row.embeddings.df[ix3.pool,1])
                          [order(apply(row.embeddings.df[ix3.pool,], 1, min_dist_to_octant_line), decreasing = TRUE)][1:n.wild])
        ix3 <- ix3[!is.na(ix3)]
        if (length(ix3)>=n.wild){
          octant3 <- row.octants[ix3]
          octant3.found <- TRUE
        } else octant3.pool <- c(octant3.pool, (octant3.pool-1) %% 8 , (octant3.pool-1) %% 8)
      }
      
      
      octant4.pool <- unique((octant2+1) %% 8 + 1)
      octant4.found <- FALSE
      
      while (!octant4.found){
        ix4.pool <- as.numeric(names(row.octants[row.octants %in% octant4.pool]))
        ix4 <- as.numeric(names(row.embeddings.df[ix4.pool,1])
                          [order(apply(row.embeddings.df[ix4.pool,], 1, min_dist_to_octant_line), decreasing = TRUE)][1:n.wild])
        ix4 <- ix4[!is.na(ix4)]
        if (length(ix4)>=n.wild){
          octant4 <- row.octants[ix4]
          octant4.found <- TRUE
        } else octant4.pool <- c(octant4.pool, (octant4.pool-1) %% 8 , (octant4.pool-1) %% 8)
      }
      
      fixed_row_index <- unique(c(ix1, ix2, ix3, ix4))
      N_fixed_row <- length(fixed_row_index)
      
      ### Find the Four Wildest Col Actors, prefarbly in Different Octants
      
      ix1 <- order(apply(col.embeddings.matrix, 1, min_dist_to_octant_line),decreasing=T)[1]
      octant1 <- col.octants[ix1]
      
      octant1.pool <- octant1
      octant1.found <- FALSE
      
      while (!octant1.found){
        ix1.pool <- as.numeric(names(col.octants[col.octants %in% octant1.pool]))
        ix1 <- as.numeric(names(col.embeddings.matrix[ix1.pool,1])
                          [order(apply(col.embeddings.matrix[ix1.pool,], 1, min_dist_to_octant_line), decreasing = TRUE)][1:n.wild])
        ix1 <- ix1[!is.na(ix1)]
        if (length(ix1)>=n.wild){
          octant1 <- col.octants[ix1]
          octant1.found <- TRUE
        } else octant1.pool <- c(octant1.pool, (octant1.pool-1) %% 8 , (octant1.pool+1) %% 8)
      }
      
      octant2.pool <- unique(octant1 + 4) 
      octant2.found <- FALSE
      
      while (!octant2.found){
        ix2.pool <- as.numeric(names(col.octants[col.octants %in% octant2.pool]))
        ix2 <- as.numeric(names(col.embeddings.matrix[ix2.pool,1])
                          [order(apply(col.embeddings.matrix[ix2.pool,], 1, min_dist_to_octant_line), decreasing = TRUE)][1:n.wild])
        ix2 <- ix2[!is.na(ix2)]
        if (length(ix2)>=n.wild){
          octant2 <- col.octants[ix2]
          octant2.found <- TRUE
        } else octant2.pool <- c(octant2.pool, (octant2.pool-1) %% 8 , (octant2.pool+1) %% 8)
      }
      
      octant3.pool <- unique((octant1+1) %% 8 + 1)
      octant3.found <- FALSE
      
      while (!octant3.found){
        ix3.pool <- as.numeric(names(col.octants[col.octants %in% octant3.pool]))
        ix3 <- as.numeric(names(col.embeddings.matrix[ix3.pool,1])
                          [order(apply(col.embeddings.matrix[ix3.pool,], 1, min_dist_to_octant_line), decreasing = TRUE)][1:n.wild])
        ix3 <- ix3[!is.na(ix3)]
        if (length(ix3)>=n.wild){
          octant3 <- col.octants[ix3]
          octant3.found <- TRUE
        } else octant3.pool <- c(octant3.pool, (octant3.pool-1) %% 8 , (octant3.pool-1) %% 8)
      }
      
      
      octant4.pool <- unique((octant2+1) %% 8 + 1)
      octant4.found <- FALSE
      
      while (!octant4.found){
        ix4.pool <- as.numeric(names(col.octants[col.octants %in% octant4.pool]))
        ix4 <- as.numeric(names(col.embeddings.matrix[ix4.pool,1])
                          [order(apply(col.embeddings.matrix[ix4.pool,], 1, min_dist_to_octant_line), decreasing = TRUE)][1:n.wild])
        ix4 <- ix4[!is.na(ix4)]
        if (length(ix4)>=n.wild){
          octant4 <- col.octants[ix4]
          octant4.found <- TRUE
        } else octant4.pool <- c(octant4.pool, (octant4.pool-1) %% 8 , (octant4.pool-1) %% 8)
      }
      
      fixed_col_index <- unique(c(ix1, ix2, ix3, ix4))
      N_fixed_col <- length(fixed_col_index)
      
      }
    
    else if(choose.method=="axis"){
      row.max <- n.wild
      col.max <- n.wild
      
      max.row.x.pool <- df_fit.mean[grep("^row_embedding\\[\\d+,1\\]$", names(df_fit.mean))]
      max.row.x.index <- as.numeric(gsub("^row_embedding\\[(\\d+),1\\]","\\1",names(max.row.x.pool[order(max.row.x.pool,decreasing=T)[1:row.max]])))
      min.row.x.index <- as.numeric(gsub("^row_embedding\\[(\\d+),1\\]","\\1",names(max.row.x.pool[order(max.row.x.pool,decreasing=F)[1:row.max]])))
      
      max.row.y.pool <- df_fit.mean[grep("^row_embedding\\[\\d+,2\\]$", names(df_fit.mean))]
      max.row.y.index <- as.numeric(gsub("^row_embedding\\[(\\d+),2\\]","\\1",names(max.row.y.pool[order(max.row.y.pool,decreasing=T)[1:row.max]])))
      min.row.y.index <- as.numeric(gsub("^row_embedding\\[(\\d+),2\\]","\\1",names(max.row.y.pool[order(max.row.y.pool,decreasing=F)[1:row.max]])))
      
      fixed_row_index <- sort(unique(c(max.row.x.index, min.row.x.index,
                                       max.row.y.index, min.row.y.index)))
      N_fixed_row <- length(fixed_row_index)
      
      max.col.x.pool <- df_fit.mean[grep("^col_embedding\\[1,\\d+\\]$", names(df_fit.mean))]
      max.col.x.index <- as.numeric(gsub("^col_embedding\\[1,(\\d+)\\]","\\1",names(max.col.x.pool[order(max.col.x.pool,decreasing=T)[1:col.max]])))
      min.col.x.index <- as.numeric(gsub("^col_embedding\\[1,(\\d+)\\]","\\1",names(max.col.x.pool[order(max.col.x.pool,decreasing=F)[1:col.max]])))
      
      max.col.y.pool <- df_fit.mean[grep("^col_embedding\\[2,\\d+\\]$", names(df_fit.mean))]
      max.col.y.index <- as.numeric(gsub("^col_embedding\\[2,(\\d+)\\]","\\1",names(max.col.y.pool[order(max.col.y.pool,decreasing=T)[1:col.max]])))
      min.col.y.index <- as.numeric(gsub("^col_embedding\\[2,(\\d+)\\]","\\1",names(max.col.y.pool[order(max.col.y.pool,decreasing=F)[1:col.max]])))
      
      fixed_col_index <- sort(unique(c(max.col.x.index, min.col.x.index,
                                       max.col.y.index, min.col.y.index)))
      N_fixed_col <- length(fixed_col_index)
      
    }
    
    fixed_row_embedding <- row.embeddings.df[fixed_row_index,]
    
    fixed_col_embedding <- col.embeddings.df[fixed_col_index,]
    fixed_col_embedding <- t(fixed_col_embedding)
  }
  
  res <- list(N_fixed_row=N_fixed_row, fixed_row_index=fixed_row_index, fixed_row_embedding=fixed_row_embedding,
              N_fixed_col=N_fixed_col, fixed_col_index=fixed_col_index, fixed_col_embedding=fixed_col_embedding)
  
  class(res) <- "LSNM_fixed_actors"
  
  return(res)
}

find.octants <- function(v){
  if (v[1] > 0 & v[2] > 0){
    if (abs(v[1]) > abs(v[2])) return (1)
    else return(2)
  }
  
  if (v[1] < 0 & v[2] > 0){
    if (abs(v[1]) > abs(v[2])) return (4)
    else return(3)
  }
  
  if (v[1] < 0 & v[2] < 0){
    if (abs(v[1]) > abs(v[2])) return (5)
    else return(6)
  } 
  
  if (v[1] > 0 & v[2] < 0){
    if (abs(v[1]) > abs(v[2])) return (8)
    else return(7)
  }
}

min_dist_to_octant_line <- function(v){
  dist.y <- abs(v[1])
  dist.x <- abs(v[2])
  dist.x.y <- abs(v[1]+v[2])/sqrt(2)
  dist.x.n.y <- abs(v[1]-v[2])/sqrt(2)
  return (min(dist.y, dist.x, dist.x.y, dist.x.n.y))
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
