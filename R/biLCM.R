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

#'@export biLCM
#'@examples \dontrun{
#'set.seed(11)
#'sim.data <- random_biLCM_data(m = 100, n = 50, k = 4, kappa_weight = NULL, a = 10000, b = 1, alpha_membership = NULL, alpha_c = rep(0.5, 100), beta_membership = NULL, beta_c = rep(0.5, 50), non_zero = TRUE)
#'res <- biLCM(edges = sim.data$A, group1.id = NULL, group2.id = NULL, count.id = NULL, k = 4, tolerance = 1e-6, max.iter = 200)
#'plot.compare.biLCM(res, sim.data, group1 = TRUE, nth = 10)
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
    stop("'edges' should be provided.\n")
  if (!class(edges)%in%c("matrix","data.frame", "igraph"))
    stop("'edges' should be matrix or data.frame or igraph.\n")
  if (class(edges)=="data.frame"&is.null(group1.id))
    stop("'group1.id' should be provided.\n")
  if (class(edges)=="data.frame"&is.null(group2.id))
    stop("'group2.id' should be provided.\n")
  if (class(edges)%in%c("data.frame", "igraph")&is.null(count.id))
    stop("'count.id' should be provided.\n")
  if (is.null(k))
    stop("'k' should be provided.\n")
  
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
    q[which(is.na(q))] <- 0
    
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

#'@name plot biLCM
#'@param biLCM_Object A trained object of class biLCM
#'@param group1 If group1 = TRUE, show the result of a member of group1 (i.e., alpha_iz's). Otherwise that of group2 (i.e., beta_jz's)
#'@param nth Show the result of 'n'th member's community distribution
#'@return A pie chart of estimated community distribution (either alpha_iz's or beta_jz's)

#'@export plot.biLCM

plot.biLCM <- function(biLCM_Object,
                       group1 = TRUE,
                       nth) {
  if(class(biLCM_Object)!="biLCM") stop("'biLCM_Object' is not of class 'biLCM'.\n")
  par(mfrow=c(1,2))
  
  k <- length(biLCM_Object$kappa)
  cols <- gg_color_hue(k, alpha = 0.7)
  
  if (group1) {
    pie(sort(biLCM_Object$alpha[nth,]), order(biLCM_Object$alpha[nth,]), col = cols, clockwise=TRUE,
        main = paste0("Estimated\nCommunity Distribution (i=",nth,")"))
  } else {
    pie(sort(biLCM_Object$beta[nth,]), order(biLCM_Object$beta[nth,]), col = cols, clockwise=TRUE,
        main = paste0("Estimated\nCommunity Distribution (j=",nth,")"))
  }
  
  par(mfrow=c(1,1))
  
}

#'Plot biLCM with latent position
#'
#'@param biLCM_Object A trained object of class biLCM
#'@param group1_space A matrix representing the latent group1 space. This matrix should have rows equal to the number of group1, and columns equal to the dimensionality of the latent space, either 1 or 2.
#'@param group2_space A matrix representing the latent group2 space. This matrix should have rows equal to the number of group2, and columns equal to the dimensionality of the latent space, either 1 or 2.
#'@param LSNM_Object A trained object of class LSNM. An optional argument, required if \code{group1_latent_position} and \code{group2_latent_position} are missing
#'@param main
#'@param legend_position
#'@param radius
#'@return plot
#'@import ggplot2
#'@import scatterpie
#'@export plot.biLCM.position

plot.biLCM.position <- function(biLCM_Object,
                                group1_latent_position = NULL,
                                group2_latent_position = NULL,
                                LSNM_Object = NULL,
                                main = "biLCM Community Distribution",
                                legend_position = "topleft",
                                radius = 0.1){
  if(class(biLCM_Object)!="biLCM") stop("'biLCM_Object' is not of class 'biLCM'.\n")
  
  m <- dim(biLCM_Object$alpha)[1] # number of group1
  n <- dim(biLCM_Object$beta)[1] # number of group2
  
  if((is.null(group1_latent_position)|is.null(group2_latent_position))&is.null(LSNM_Object)) stop("Either group1_latent_position/group2_latent_position or LSNM_Object should be provided.\n")
  
  if(is.null(LSNM_Object)) {
    if(dim(group1_latent_position)[1]!=m) stop("Invalid number of rows in 'group1_latent_position'.\n")
    if(dim(group2_latent_position)[1]!=n) stop("Invalid number of rows in 'group2_latent_position'.\n")
    
    D <- dim(group1_latent_position)[2]
    
    if(dim(group2_latent_position)[2]!=D) stop("Invalid number of cols in either 'group1_latent_position' or 'group2_latent_position'.\n")
    
    if(D==1) {
      dat <- data.frame(x = c(group1_latent_position[,1],
                              group2_latent_position[,1]),
                        name = c(paste0("group1.",1:m), paste0("group2.",1:n)),
                        group = c(rep("group1",m), rep("group2",n)))
      dat$y <- match(dat$x, sort(dat$x, decreasing = T)) 
      dat$y <- (dat$y/max(dat$y))*(range(dat$x)[2]-range(dat$x)[1])
      
      community <- as.data.frame(rbind(biLCM_Object$alpha, biLCM_Object$beta))
      colnames(community) <- LETTERS[1:length(biLCM_Object$kappa)]
      dat <- cbind(dat, community)
      
      ggplot() + geom_scatterpie(aes(x=x, y=y, group=group, r=radius), data=dat,
                                 cols=LETTERS[1:length(biLCM_Object$kappa)], color=NA, alpha = 0.7) + 
        coord_equal() + theme_classic() + theme(legend.position = legend_position, plot.title = element_text(hjust=0.5)) + 
        scale_y_continuous(breaks=dat$y, labels=dat$name) +
        xlab("Latent Space Dimension 1") + ylab(NULL) +
        labs(title = main, fill = "Community")
    } else if (D==2) {
      dat <- data.frame(x = c(group1_latent_position[,1],
                              group2_latent_position[,1]),
                        y = c(group1_latent_position[,2],
                              group2_latent_position[,2]),
                        group = c(rep("group1",m), rep("group2",n)))
      community <- as.data.frame(rbind(biLCM_Object$alpha, biLCM_Object$beta))
      colnames(community) <- LETTERS[1:length(biLCM_Object$kappa)]
      dat <- cbind(dat, community)
      
      ggplot() + geom_scatterpie(aes(x=x, y=y, group=group, r=radius), data=dat,
                                 cols=LETTERS[1:length(biLCM_Object$kappa)], color=NA, alpha = 0.7) + 
        coord_equal() + theme_classic() + theme(legend.position = legend_position, plot.title = element_text(hjust=0.5)) + 
        xlab("Latend Space Dimension 1") + ylab("Latend Space Dimension 2") +
        labs(title = main, fill = "Community")
    }
    
  } else {
    if(LSNM_Object$stan_fitted_model@par_dims$row_factor_adj!=m) stop("Invalid number of group1 members in 'LSNM_Object'")
    if(LSNM_Object$stan_fitted_model@par_dims$col_factor_adj!=n) stop("Invalid number of group2 members in 'LSNM_Object'")
    
    D <- LSNM_Object$stan_fitted_model@par_dims$cov_embedding_diag # number of dimensions
    
    df_fit <- as.data.frame(LSNM_Object$stan_fitted_model)
    nms <- df_fit[ , grepl( "^col_embedding|^row_embedding|^col_factor|^row_factor" , names(df_fit) )]
    plot.data <- colMeans(nms) # posterior mean
    
    if (D==1) {
      dat <- data.frame(x = c(plot.data[paste0("row_embedding[",1:m,",1]")],
                              plot.data[paste0("col_embedding[1,",1:n,"]")]),
                        name = c(paste0("group1.",1:m), paste0("group2.",1:n)),
                        group = c(rep("group1",m), rep("group2",n)))
      dat$y <- match(dat$x, sort(dat$x, decreasing = T)) 
      dat$y <- (dat$y/max(dat$y))*(range(dat$x)[2]-range(dat$x)[1])
      
      community <- as.data.frame(rbind(biLCM_Object$alpha, biLCM_Object$beta))
      colnames(community) <- LETTERS[1:length(biLCM_Object$kappa)]
      dat <- cbind(dat, community)
      
      ggplot() + geom_scatterpie(aes(x=x, y=y, group=group, r=radius), data=dat,
                                 cols=LETTERS[1:length(biLCM_Object$kappa)], color=NA, alpha = 0.7) + 
        coord_equal() + theme_classic() + theme(legend.position = legend_position, plot.title = element_text(hjust=0.5)) + 
        scale_y_continuous(breaks=dat$y, labels=dat$name) +
        xlab("LSNM Dimension 1") + ylab(NULL) +
        labs(title = main, fill = "Community")
      
    } else if (D==2) {
      dat <- data.frame(x = c(plot.data[paste0("row_embedding[",1:m,",1]")],
                              plot.data[paste0("col_embedding[1,",1:n,"]")]),
                        y = c(plot.data[paste0("row_embedding[",1:m,",2]")],
                              plot.data[paste0("col_embedding[2,",1:n,"]")]),
                        group = c(rep("group1",m), rep("group2",n)))
      community <- as.data.frame(rbind(biLCM_Object$alpha, biLCM_Object$beta))
      colnames(community) <- LETTERS[1:length(biLCM_Object$kappa)]
      dat <- cbind(dat, community)
      
      ggplot() + geom_scatterpie(aes(x=x, y=y, group=group, r=radius), data=dat,
                                 cols=LETTERS[1:length(biLCM_Object$kappa)], color=NA, alpha = 0.7) + 
        coord_equal() + theme_classic() + theme(legend.position = legend_position, plot.title = element_text(hjust=0.5)) + 
        xlab("LSNM Dimension 1") + ylab("LSNM Dimension 2") +
        labs(title = main, fill = "Community")
    }
  }
}
