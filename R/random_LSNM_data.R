#'@param D The dimensionality of the latent space, 2 dimensions is recommended
#'@param group1 The number of group1
#'@param group2 The number of group2
#'@param group1_space The desired latent space positions of the group1. Expecting a matrix with group1 rows and D columns. Assumed to be drawn from a multivariate normal distribution with a mean of 0. An optional argument, but must specify this or tau.
#'@param group2_space The desired latent space positions of the group2. Expecting a matrix with group2 rows and D columns. Assumed to be drawn from a multivariate normal distribution. An optional argument, but must specifiy this or mu and Sigma.
#'@param mu The mean vector to generate group2_space. Expecting a vector of length D. An optional argument, required if missing group2_space.
#'@param Sigma The variance-covariance matrix used to generate group2_space. Expecting a D by D matrix with a reasonable variance-covariance structure. An optional argument, required if missing group2_space.
#'@param tau The variance-covariance matrix used to generate group1_space. Expecting a D by D matrix with a diagonal variance-covariance structure. An optional argument, required if missing group1_space.
#'@param alpha_popularity The group1 popularity factors used to account for baseline likelihood to lobby. Expecting a group1 length vector. Assumed to be drawn from a normal distribution. An optional argument, but must specify this or v and sigma_sq_L
#'@param Beta_popularity The group2 popularity factors used to account for baseline likelihood to sponsor bills. Expecting a group2 length vector. Assumed to be drawn from a normal distribution with mean 0. An optional argument, but must specify this or sigma_sq_P
#'@param v The mean used to generate the alpha_popularity vector. An optional argument, required if missing alpha_popularity
#'@param sigma_sq_L The variance used to generate the alpha_popularity vector. An optional argument, required if missing alpha_popularity
#'@param sigma_sq_P The variance used to generate the Beta_popularity vector. An optional argument, required if missing Beta_popularity
#'@return A list of the latent group1 space, the latent group2 space, and then the randomly generated poisson count matrix A.

#'@import MASS
#'@useDynLib polnet, .registration = TRUE
#'@export random_LSNM_data


random_LSNM_data <- function(D,
                             group1,
                             group2,
                             group1_space = NULL,
                             group2_space = NULL,
                             mu = NULL,
                             Sigma = NULL,
                             tau = NULL,
                             alpha_popularity = NULL,
                             Beta_popularity = NULL,
                             v = NULL,
                             sigma_sq_L = NULL,
                             sigma_sq_P = NULL){
  #Generate the alpha, Beta vectors if not existing already
  if(is.null(alpha_popularity)){
    if(is.null(v) | is.null(sigma_sq_L)){
      stop("Invalid Alpha Parameters")
    }else{
      alpha_popularity = rnorm(group1, v, sigma_sq_L)
    }
  }
  if(is.null(Beta_popularity)){
    if(is.null(sigma_sq_P)){
      stop("Invalid Beta Parameters")
    }else{
      Beta_popularity = rnorm(group2, 0, sigma_sq_P)
    }
  }
  #putting together matrix of alpha + Beta
  a_B = outer(alpha_popularity, Beta_popularity, '+')
  if(is.null(group1_space)){
    if(is.null(tau)){
      stop('Invalid Theta Parameters')
    }else{
      if(!all(tau[lower.tri(tau)] == 0, tau[upper.tri(tau)] == 0)){
        stop('Tau covariance matrixmust be diagonal')
      }
      group1_space = mvrnorm(n = group1, mu = rep(0, D), Sigma = tau)
    }
  }
  if(is.null(group2_space)){
    if(is.null(mu) | is.null(Sigma)){
      stop('Invalid Psi Parameters')
    }else{
      group2_space = mvrnorm(n = group2, mu = mu, Sigma = Sigma)
    }
  }
  if (D>=2) distances = outer(split(group1_space, row(group1_space)), split(group2_space, row(group2_space)), Vectorize(l2_norm.sq))
	else distances = outer(group1_space, group2_space, Vectorize(l2_norm.sq))
  diff_vec = a_B - distances
  A_mat = matrix(sapply(as.vector(diff_vec), function(x) rpois(1, exp(x))), nrow = nrow(distances))
  return(list(Theta = group1_space, Psi = group2_space, A = A_mat))
}


l1_norm <- function(vec1, vec2){
  return(sum(abs(vec1 - vec2)))
}

l2_norm <- function(vec1, vec2){
  return(sqrt(sum((vec1 - vec2)^2)))
}

l2_norm.sq <- function(vec1, vec2){
  return(sum((vec1 - vec2)^2))
}

#'@param m
#'@param n
#'@param k
#'@param m.center.cor
#'@param n.center.cor
#'@return A list of simulated data

#'@import MASS
#'@useDynLib polnet, .registration = TRUE
#'@export

# Generate clustered LSNM data
random_LSNM_data_cluster <- function(n.cluster=4,
                                     group1.center=rbind(c(-1,-1), c(-1, 1), c(1, 1), c(1, -1))*5,
                                     group2.center=rbind(c(-1,-1), c(-1, 1), c(1, 1), c(1, -1))*5,
                                     D=2,
                                     group1=100,
                                     group2=50,
                                     Sigma = NULL,
                                     tau = 1,
                                     v = 0,
                                     sigma_sq_L = 2,
                                     sigma_sq_P = 3){

  # Generate Positive Definite Matrices
  Posdef <- function (n, ev = runif(n, 0, 1)) {
    Z <- matrix(ncol=n, rnorm(n^2))
    decomp <- qr(Z)
    Q <- qr.Q(decomp)
    R <- qr.R(decomp)
    d <- diag(R)
    ph <- d / abs(d)
    O <- Q %*% diag(ph)
    Z <- t(O) %*% diag(ev) %*% O
    return(Z)
  }

  if (is.null(Sigma)){
    if (D>=2){
      Sigma <- Posdef(n=D)
    }
    else Sigma <- runif(1, 1, 2)
  }

  # Initialize Actor Coordinates
  group1.cor <- matrix(0, nrow=0, ncol=D+1)
  group2.cor <- matrix(0, nrow=0, ncol=D+1)

  # Generate number of units in each clusters randomly
  group1.division <- rmultinom(n=1, size=group1, prob=rep(1/n.cluster,n.cluster))
  group2.division <- rmultinom(n=1, size=group2, prob=rep(1/n.cluster,n.cluster))

  # Simulated Coordinates Generation
  if (n.cluster>=2){
    for (i in 1:n.cluster){
      # Generate Coordinates for Actor 1
      if (D>=2){
        group1.sigma <- diag(D)*tau
      } else group1.sigma <- tau
      temp.new.group1.cor <- mvrnorm(n=group1.division[i], mu=group1.center[i,],
                                      Sigma=group1.sigma)

      temp.new.group1.cor <- cbind(temp.new.group1.cor, i)
      group1.cor <- rbind(group1.cor, temp.new.group1.cor)

      # Generate Coordinates for Actor 2
      group2.sigma <- Sigma
      temp.new.group2.cor <- mvrnorm(n=group2.division[i], mu=group2.center[i,],
                                          Sigma=group2.sigma)
      temp.new.group2.cor <- cbind(temp.new.group2.cor, i)
      group2.cor <- rbind(group2.cor, temp.new.group2.cor)
    }
  } else {

    # Generate Coordinates for Actor 1
    if (D>=2){
      group1.sigma <- diag(D)*tau
    } else group1.sigma <- tau
    temp.new.group1.cor <- mvrnorm(n=group1.division, mu=group1.center,
                                    Sigma=group1.sigma)

    temp.new.group1.cor <- cbind(temp.new.group1.cor, 1)
    group1.cor <- rbind(group1.cor, temp.new.group1.cor)

    # Generate Coordinates for Actor 2
    group2.sigma <- Sigma
    temp.new.group2.cor <- mvrnorm(n=group2.division, mu=group2.center,
                                        Sigma=group2.sigma)
    temp.new.group2.cor <- cbind(temp.new.group2.cor, 1)
    group2.cor <- rbind(group2.cor, temp.new.group2.cor)

  }

  # Randomly Generate Actors Popularity
  group1.popularity <- rnorm(group1, v, sqrt(sigma_sq_L))
  group2.popularity <- rnorm(group2, 0, sqrt(sigma_sq_P))


  # Call existing function
  LSNM_data <- random_LSNM_data(D=D, group1=group1, group2 = group2, group1_space = group1.cor[,-(D+1)],
                                group2_space = group2.cor[,-(D+1)], alpha_popularity = group1.popularity,
                                Beta_popularity = group2.popularity)

  # return
  ret.list <- list(LSNM_data=LSNM_data, group1.popularity=group1.popularity, group2.popularity=group2.popularity, group1.cor=group1.cor,
                   group2.cor=group2.cor, D=D, group1.center=group1.center,
                   group2.center=group2.center, group1=group1,
                   group2=group2, Sigma=Sigma, tau=tau, v=v, sigma_sq_L = sigma_sq_L,
                   sigma_sq_P = sigma_sq_P,
                   group1.cluster = rep(1:n.cluster, group1.division),
                   group2.cluster = rep(1:n.cluster, group2.division))
}

#'@param LSNM_Object A trained object of class LSNM
#'@param group1_space A matrix representing the true latent group1 space. This matrix should have rows equal to the number of group1, and columns equal to the dimensionality of the latent space, D.
#'@param group2_space A matrix representing the true latent group2 space. This matrix should have rows equal to the number of group2, and columns equal to the dimensionality of the latent space, D.
#'@return Does not return an object. Prints the proportion of latent space estimates that fell within the credible interval as well as the average error from the true latent space estimates.
#'

#'@useDynLib polnet, .registration = TRUE
#'@export


compare.LSNM <- function(LSNM_Object,
                         group1_space,
                         group2_space){
  lsnmobj <- summary.LSNM(LSNM_Object)
  l_ordered <- lsnmobj[order(rownames(lsnmobj)), ]
  tru_pars <- c(as.vector(t(group1_space)), as.vector(t(group2_space)))
  perc_in_cred <- sum(tru_pars < l_ordered$`90%` & tru_pars > l_ordered$`10%`) / length(tru_pars)
  ave_marg_error <- mean(abs(l_ordered$Mean - tru_pars))
  paste0('A proportion of ', perc_in_cred, 'of latent space parameters fell within their credible interval for an average error of ', ave_marg_error)
}

#'@param LSNM_Object A trained object of class LSNM
#'@param group1_space A matrix representing the true latent group1 space. This matrix should have rows equal to the number of group1, and columns equal to the dimensionality of the latent space, D.
#'@param group2_space A matrix representing the true latent group2 space. This matrix should have rows equal to the number of group2, and columns equal to the dimensionality of the latent space, D.
#'@param group1_popularity The group1 popularity factors used to account for baseline likelihood to lobby. Expecting a group1 length vector. Assumed to be drawn from a normal distribution. An optional argument, but must specify this or v and sigma_sq_L
#'@param group2_popularity The group2 popularity factors used to account for baseline likelihood to sponsor bills. Expecting a group2 length vector. Assumed to be drawn from a normal distribution with mean 0. An optional argument, but must specify this or sigma_sq_P
#'@param group1_cluster A vector representing the true cluster of group1
#'@param group2_cluster A vector representing the true cluster of group2
#'@return Two plots: true latent space and estimated LSNM positions

#'@useDynLib polnet, .registration = TRUE
#'@export


plot.compare.LSNM <- function(LSNM_Object,
                              group1_space,
                              group2_space,
                              group1_popularity,
                              group2_popularity,
                              group1_cluster = NULL,
                              group2_cluster = NULL,
                              main = "Estimated LSNM Positions",
                              legend = c("Group1", "Group2"),
                              legend_position = "topleft",
                              ...){
  D <- ifelse(is.null(ncol(group1_space)),1,2) # number of dimensions
  m <- length(group1_space)/D
  n <- length(group2_space)/D

  group1_col <- group1_cluster
  group2_col <- group2_cluster

  if (is.null(group1_col)) group1_col <- rep("black", m)
  if (is.null(group2_col)) group2_col <- rep("black", n)

  if (D==1) {
    df_fit <- as.data.frame(LSNM_Object$stan_fitted_model)
    nms <- df_fit[ , grepl( "^col_embedding|^row_embedding|^col_factor|^row_factor" , names(df_fit) )]
    plot.data <- colMeans(nms) # posterior mean

    row_elements <- plot.data[paste0("row_embedding[",1:m,",1]")]
    col_elements <- plot.data[paste0("col_embedding[",1:n,",1]")]

    plot(x = row_elements,
         y = group1_space,
         pch = 1,
         cex = 1,
         xlab = "Estimate Dimension 1",
         ylab = "True Dimension 1",
         yaxt = "n",
         col = group1_col,...)
    points(x = col_elements,
           y = group2_space,
           pch = 0,
           cex = 1,
           col = group2_col)

  } else {
    row_size <- exp(group1_popularity) # size of group1
    row_size <- 2*row_size/max(row_size)
    col_size <- exp(group2_popularity) # size of group2
    col_size <- 2*col_size/max(col_size)

    par(mfrow=c(1,2))

    plot(x = group1_space[,1],
         y = group1_space[,2],
         pch = 1,
         cex = row_size,
         xlab = "Latent Space Dimension 1",
         ylab = "Latent Space Dimension 2",
         main = "True Latent Space",
         col = group1_col, ...)
    points(x = group2_space[,1],
           y = group2_space[,2],
           pch = 0,
           cex = col_size,
           col = group2_col)

    plot.LSNM(LSNM_Object, group1_cluster, group2_cluster, main, legend, legend_position, ...)
    
    par(mfrow=c(1,1))
  }

}
