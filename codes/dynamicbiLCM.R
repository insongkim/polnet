library(gtools) # rdirichlet
library(mvtnorm) # rmvnorm

gg_color_hue <- function(n, alpha = 1) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100, alpha)[1:n]
}

# from https://gist.github.com/aufrank/83572
logsumexp <- function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}

softmax <- function (x) {
  exp(x - logsumexp(x))
}

# 1. Synthetic Data -------------------------------------------------------

random.dynamic.biLCM <- function(## Data size
                                 # number of community = length of z (+int)
                                 k = 3,
                                 # number of member state = length of i (+int)
                                 m = 5,
                                 # number of time period = length of t (+int; Assmp: \underbar{T_i}=1 and \upperbar{T_i}=s \forall i)
                                 s = 2,
                                 # number of bills for each time period = length of j's (+int; Assmp: n_1 = ... = n_s)
                                 n = 10,
                          
                                 ## Hyper-parameters
                                 #-- a: shape parameter to generate kappa (+real)
                                 #-- b: rate parameter to generate kappa (+real)
                                 a = 2,
                                 b = 0.5,
                                 #-- state-space model: a_it | a_i,t-1 ~ N(a_i,t-1, sigma^2 I)
                                 sigma = 2,
                                 #-- beta_c: concentration vector to generate beta membership (vector of legnth m; +real)
                                 beta_c = rep(0.5, n),
                                 
                                 ## random seed
                                 true_seed = 1
                                 ) {
  
  set.seed(true_seed)
  ## Parameters
  # kappa: the overall level of activity in community z at time t (kxs matrix; +real)
  # #-- a: shape parameter to generate kappa (+real)
  # #-- b: rate parameter to generate kappa (+real)
  # a = 2
  # b = 0.5
  kappa = matrix(rgamma(k*s, shape = a, rate = b), k, s)
  # softmax(alpha): member state i's involvement in communities at time t (mxsxk array; +real)
  #-- state-space model: a_it | a_i,t-1 ~ N(a_i,t-1, sigma^2 I)
  alpha = array(NA, dim = c(m,s,k))
  # sigma = 2
  alpha[,1,] = matrix(rnorm(m*k, 0, sigma), m, k)
  for (t in 2:s) {
    for (i in 1:m) {
      alpha[i,t,] = rmvnorm(1, mean = alpha[i,t-1,], sigma = diag(sigma^2, k))
    }
  }
  softmax_alpha = apply(alpha, c(2,3), softmax)
  # beta: bill j's involvement in communities at time t (nxsxk array; +real)
  # #-- beta_c: concentration vector to generate beta membership (vector of legnth m; +real)
  # beta_c = rep(0.5, n)
  beta = array(NA, dim = c(n,s,k))
  for (t in 1:s) {
    beta[,t,] = t(rdirichlet(k, beta_c))
  }
  
  ## Data
  # adjacency matrix (mxnxs array; 0/+int)
  A = array(NA, dim = c(m,n,s))
  for (i in 1:m) {
    for (j in 1:n) {
      for (t in 1:s) {
        A[i,j,t] = rpois(1, sum(kappa[,t] * softmax_alpha[i,t,] * beta[j,t,]))
      }
    }
  }
  
  
  ## save
  out = list(kappa = kappa,
             alpha = alpha,
             softmax_alpha = softmax_alpha,
             beta = beta,
             A = A)
  
  return(out)
}


# 2. Estimation -----------------------------------------------------------

dynamic.biLCM <- function(## Data
                          A,
                          # number of community = length of z (+int)
                          k = 3,
                          
                          ## Hyper-parameters/Initialization
                          # kappa
                          a = 2,
                          b = 0.5,
                          # alpha
                          sigma = 2,
                          # beta
                          beta_c = rep(1, n),
                          # q
                          q_c = rep(1, k),
                          # Variational parameters for updating alpha
                          # nu_hat (+real)
                          nu_hat = 2,
                          # Sigma = sigma^2*I
                          # Nu = nu_hat^2*I
                          
                          ## Tuning parameters
                          # tolerance
                          tolerance = 1e-6,
                          # step size for updating alpha_hat
                          r = 0.001,
                          # maximum iterations
                          max_iter = 5000,
                          
                          ## random seed for initialization
                          init_seed = 123,
                          
                          ## initial conditions for forward (alpha_hat)
                          # fixed mi0 = 0 and Vi0 = I
                          kappa = NULL,
                          alpha = NULL,
                          beta = NULL
                          ) {
  ## Data size
  # # number of community = length of z (+int)
  # k = 3
  # number of member state = length of i (+int)
  m = dim(A)[1]
  # # number of time period = length of t (+int; Assmp: \underbar{T_i}=1 and \upperbar{T_i}=s \forall i)
  s = dim(A)[3]
  # number of bills for each time period = length of j's (+int; Assmp: n_1 = ... = n_s)
  n = dim(A)[2]
  # number of iterations
  n_iter = 0
  # log-likelihood
  llik = rep(NA, n_iter)
  # a diagonal matrix
  delta = diag(1, s, s)
  I = diag(1, k, k)
  
  # Initialization ----------------------------------------------------------
  set.seed(init_seed)
  # # kappa
  # a = 2
  # b = 0.5
  if (is.null(kappa)) {
    kappa = matrix(rgamma(k*s, shape = a, rate = b), k, s)
  }
  # alpha
  if (is.null(alpha)) {
    alpha = array(NA, dim = c(m,s,k))
    # sigma = 2
    alpha[,1,] = matrix(rnorm(m*k, 0, sigma), m, k)
    for (t in 2:s) {
      for (i in 1:m) {
        alpha[i,t,] = rmvnorm(1, mean = alpha[i,t-1,], sigma = diag(sigma^2, k))
      }
    }
  }
  softmax_alpha = apply(alpha, c(2,3), softmax)
  # beta
  # beta_c = rep(1, n)
  if (is.null(beta)) {
    beta = array(NA, dim = c(n,s,k))
    for (t in 1:s) {
      beta[,t,] = t(rdirichlet(k, beta_c))
    }
  }
  # q (mxnxsxk array)
  # q = array(0, c(m,n,s,k))
  # q_c = rep(1, k)
  q = array(NA, dim = c(m,n,s,k))
  for (i in 1:m) {
    for (j in 1:n) {
      for (t in 1:s) {
        q[i,j,t,] = rdirichlet(1, q_c)
      }
    }
  }
  
  # Compute log-likelihood
  kab = array(NA, c(m,n,s))
  for (t in 1:s) {
    kab[,,t] = sweep(softmax_alpha[,t,], 2, kappa[,t], "*") %*% t(beta[,t,])
  }
  llik_old = sum(A * log(kab) - kab)
  
  # Variational parameters for updating alpha
  # nu_hat (+real)
  # nu_hat = 2
  # alpha_hat (mxsxk array)
  alpha_hat = array(NA, dim = c(m,s,k))
  for (t in 1:s) {
    for (i in 1:m) {
      alpha_hat[i,t,] = rmvnorm(1, mean = alpha[i,t,], sigma = diag(nu_hat^2, k))
    }
  }
  # compute m, V, \tilde m, \tilde V (forward/backward)
  Sigma = sigma^2*I
  Nu = nu_hat^2*I
  #-- forward
  M = array(NA, dim = c(m,s,k))
  V = array(NA, dim = c(m,s,k,k))
  # initial condition: fixed mi0 and Vi0
  M0 = matrix(0, m, k) 
  V0 = array(NA, dim = c(m,k,k)) 
  for (i in 1:m) {
    V0[i,,] = diag(1,k,k)
  }
  for (i in 1:m) {
    P1 = V0[i,,]+Sigma
    P2 = P1%*%solve(P1+Nu)
    P3 = I-P2
    M[i,1,] = P3%*%M0[i,] + P2%*%alpha_hat[i,1,]
    V[i,1,,] = P3%*%P1
  }
  for (t in 2:s) {
    for (i in 1:m) {
      P1 = V[i,t-1,,]+Sigma
      P2 = P1%*%solve(P1+Nu)
      P3 = I-P2
      M[i,t,] = P3%*%M[i,t-1,] + P2%*%alpha_hat[i,t,]
      V[i,t,,] = P3%*%P1
    }
  }
  #-- backward
  M_tilde = array(NA, dim = c(m,s,k))
  V_tilde = array(NA, dim = c(m,s,k,k))
  M_tilde_lag = array(NA, dim = c(m,s,k)) # for step 4-1
  # initial condition
  M_tilde[,s,] = M[,s,]
  V_tilde[,s,,] = V[,s,,]
  for (i in 1:m) {
    for (t in (s-1):1) {
      P1 = V[i,t,,] + Sigma
      P2 = V[i,t,,]%*%solve(P1)
      M_tilde[i,t,] = (I-P2)%*%M[i,t,] + P2%*%M_tilde[i,t+1,]
      V_tilde[i,t,,] = V[i,t,,] + P2%*%(V_tilde[i,t+1,,]-P1)%*%t(P2)
      M_tilde_lag[i,t+1,] = M_tilde[i,t,]
    }
    P1 = V0[i,,]+Sigma
    P2 = V0[i,,]%*%solve(P1)
    M_tilde_lag[i,1,] = (I-P2)%*%M0[i,] + P2%*%M_tilde[i,1,]
  }
  # zeta_hat (sxk matrix)
  zeta_hat = matrix(NA, s, k)
  for (t in 1:s) {
    for (z in 1:k) {
      zeta_hat[t,z] = sum(exp(M_tilde[,t,z]+V_tilde[,t,z,z]/2))
    }
  }
  # compute \partial l.b./\partial alpha_hat
  der = array(NA, dim = c(m, s, k))
  #-- forward: \partial m/\partial alpha_hat
  der_for = array(NA, dim = c(m, s, s, k, k))
  # initial condition
  for (i in 1:m) {
    for (t2 in 1:s) {
      P1 = V0[i,,]+Sigma
      P2 = P1%*%solve(P1+Nu)
      der_for[i,1,t2,,] = P2*delta[t2,1]
    }
  }
  for (i in 1:m) {
    for (t1 in 2:s) {
      for (t2 in 1:s) {
        P1 = V[i,t-1,,]+Sigma
        P2 = P1%*%solve(P1+Nu)
        P3 = I-P2
        der_for[i,t1,t2,,] = P3%*%der_for[i,t1-1,t2,,] + P2*delta[t2,t1]
      }
    }
  }
  #-- backward: \partial m_tilde/\partial alpha_hat
  der_back = array(NA, dim = c(m, s, s, k, k))
  der_back_lag = array(NA, dim = c(m, s, s, k, k))
  # initial condition
  for (i in 1:m) {
    for (t2 in 1:s) {
      der_back[i,s,t2,,] = der_for[i,s,t2,,]
    }
  }
  for (i in 1:m) {
    for (t1 in (s-1):1) {
      for (t2 in 1:s) {
        P1 = V[i,t1,,] + Sigma
        P2 = V[i,t1,,]%*%solve(P1)
        der_back[i,t1,t2,,] = (I-P2)%*%der_for[i,t1,t2,,] + P2%*%der_back[i,t1+1,t2,,]
        der_back_lag[i,t1+1,t2,,] = der_back[i,t1,t2,,]
      }
      for (t2 in 1:s) {
        P1 = V0[i,,]+Sigma
        P2 = V0[i,,]%*%solve(P1)
        der_back_lag[i,1,t2,,] = P2%*%der_back[i,1,t2,,]
      }
    }
  }
  # compute derivative
  for (i in 1:m) {
    for (t2 in 1:s) {
      for (z in 1:k) {
        der[i,t2,z] = -(1/sigma^2)*(M_tilde[i,,z] - M_tilde_lag[i,,z])%*%(der_back[i,,t2,z,z] - der_back_lag[i,,t2,z,z]) -
          colSums(A[i,,]*q[i,,,z] - t((t(apply(A[,,]*q[,,,z], c(2,3), sum)) * 1/zeta_hat[,z]) * colSums(exp(M_tilde[,,z] + V_tilde[,,z,z]/2)))) %*% der_back[i,,t2,z,z]
      }
    }
  }
  
  repeat {
    n_iter = n_iter + 1
    # Step 1: Update q --------------------------------------------------------
    for (t in 1:s) {
      log_numer = log(sapply(1:k, function(z) kappa[z,t]*outer(softmax_alpha[,t,z], beta[,t,z], "*"), simplify = "array"))
      log_numer[which(is.infinite(log_numer))] <- -1.797693e+308
      log_denom = apply(log_numer, c(1,2), logsumexp)
      log_denom[which(is.infinite(log_denom))] <- -1.797693e+308
      q[,,t,] = exp(sweep(log_numer, c(1,2), log_denom, "-"))
    }
    
    # Step 2: Update kappa ----------------------------------------------------
    for (t in 1:s) {
      numer <- sapply(1:k, function(z) sum(A[,,t] * q[,,t,z]))
      kappa[,t] <- numer
      # denom. always sums to 1: sapply(1:k, function(z) sum(outer(softmax_alpha[,t,z], beta[,t,z], "*")))
    }
    
    # Step 3: Update beta -----------------------------------------------------
    for (t in 1:s) {
      log_numer <- log(sapply(1:k, function(z) colSums(A[,,t] * q[,,t,z])))
      log_numer[which(is.infinite(log_numer))] <- -1.797693e+308
      log_denom <- log(sapply(1:k, function(z) sum(A[,,t] * q[,,t,z])))
      log_denom[which(is.infinite(log_denom))] <- -1.797693e+308
      beta[,t,] <- exp(sweep(log_numer, 2, log_denom, "-"))
    }
    
    # Step 4: Update alpha ----------------------------------------------------
    # 4-1: Compute m, V, \tilde m, \tilde V (forward/backward)
    #-- forward
    for (i in 1:m) {
      P1 = V0[i,,]+Sigma
      P2 = P1%*%solve(P1+Nu)
      P3 = I-P2
      M[i,1,] = P3%*%M0[i,] + P2%*%alpha_hat[i,1,]
      V[i,1,,] = P3%*%P1
    }
    for (t in 2:s) {
      for (i in 1:m) {
        P1 = V[i,t-1,,]+Sigma
        P2 = P1%*%solve(P1+Nu)
        P3 = I-P2
        M[i,t,] = P3%*%M[i,t-1,] + P2%*%alpha_hat[i,t,]
        V[i,t,,] = P3%*%P1
      }
    }
    #-- backward
    M_tilde[,s,] = M[,s,]
    V_tilde[,s,,] = V[,s,,]
    for (i in 1:m) {
      for (t in (s-1):1) {
        P1 = V[i,t,,] + Sigma
        P2 = V[i,t,,]%*%solve(P1)
        M_tilde[i,t,] = (I-P2)%*%M[i,t,] + P2%*%M_tilde[i,t+1,]
        V_tilde[i,t,,] = V[i,t,,] + P2%*%(V_tilde[i,t+1,,]-P1)%*%t(P2)
        M_tilde_lag[i,t+1,] = M_tilde[i,t,]
      }
      P1 = V0[i,,]+Sigma
      P2 = V0[i,,]%*%solve(P1)
      M_tilde_lag[i,1,] = (I-P2)%*%M0[i,] + P2%*%M_tilde[i,1,]
    }
    
    # 4-2: Update variational parameter zeta_hat
    for (t in 1:s) {
      for (z in 1:k) {
        zeta_hat[t,z] = sum(exp(M_tilde[,t,z]+V_tilde[,t,z,z]/2))
      }
    }
    
    # 4-3: Update variational parameter alpha_hat (forward/backward -> compute derivative -> old + stepsize*derivative)
    #-- forward: \partial m/\partial alpha_hat
    for (i in 1:m) {
      for (t1 in 2:s) {
        for (t2 in 1:s) {
          P1 = V[i,t-1,,]+Sigma
          P2 = P1%*%solve(P1+Nu)
          P3 = I-P2
          der_for[i,t1,t2,,] = P3%*%der_for[i,t1-1,t2,,] + P2*delta[t2,t1]
        }
      }
    }
    #-- backward: \partial m_tilde/\partial alpha_hat
    for (i in 1:m) {
      for (t1 in (s-1):1) {
        for (t2 in 1:s) {
          P1 = V[i,t1,,] + Sigma
          P2 = V[i,t1,,]%*%solve(P1)
          der_back[i,t1,t2,,] = (I-P2)%*%der_for[i,t1,t2,,] + P2%*%der_back[i,t1+1,t2,,]
          der_back_lag[i,t1+1,t2,,] = der_back[i,t1,t2,,]
        }
        for (t2 in 1:s) {
          P1 = V0[i,,]+Sigma
          P2 = V0[i,,]%*%solve(P1)
          der_back_lag[i,1,t2,,] = P2%*%der_back[i,1,t2,,]
        }
      }
    }
    # compute derivative
    for (i in 1:m) {
      for (t2 in 1:s) {
        for (z in 1:k) {
          der[i,t2,z] = -(1/sigma^2)*(M_tilde[i,,z] - M_tilde_lag[i,,z])%*%(der_back[i,,t2,z,z] - der_back_lag[i,,t2,z,z]) -
            colSums(A[i,,]*q[i,,,z] - t((t(apply(A[,,]*q[,,,z], c(2,3), sum)) * 1/zeta_hat[,z]) * colSums(exp(M_tilde[,,z] + V_tilde[,,z,z]/2)))) %*% der_back[i,,t2,z,z]
        }
      }
    }
    # old + stepsize*derivative
    alpha_hat = alpha_hat + r*der
    
    # 4-4: alpha = \tilde m
    alpha = M_tilde
    softmax_alpha = apply(alpha, c(2,3), softmax)
    
    # Compute log-likelihood --------------------------------------------------
    for (t in 1:s) {
      kab[,,t] = sweep(softmax_alpha[,t,], 2, kappa[,t], "*") %*% t(beta[,t,])
    }
    log_kab <- log(kab)
    log_kab[which(is.na(log_kab))] <- 0
    log_kab[which(is.infinite(log_kab))] <- -1.797693e+308
    llik_new <- sum(A * log_kab - kab)
    llik[n_iter] <- llik_new
    
    
    
    # Check convergence -------------------------------------------------------
    if (llik_new - llik_old < tolerance || n_iter == max_iter) break
    llik_old <- llik_new
  }
  out <- list(loglikelihood = llik, kappa = kappa, alpha = alpha, softmax_alpha = softmax_alpha, beta = beta)
}




# 3. Compare estimates ----------------------------------------------------

plot.compare.dynbiLCM <- function(dynbiLCM_Object,
                                 simulated_data,
                                 group1 = TRUE,
                                 i,
                                 t) {
  # if(class(biLCM_Object)!="biLCM") stop("'biLCM_object' is not of class 'biLCM'.\n")
  
  par(mfrow=c(1,2))
  
  k1 <- length(simulated_data$kappa)
  cols1 <- gg_color_hue(k1, alpha = 0.7)
  
  k2 <- length(dynbiLCM_Object$kappa)
  cols2 <- gg_color_hue(k2, alpha = 0.7)
  
  if (group1) {
    pie(sort(simulated_data$softmax_alpha[i,t,]), order(simulated_data$softmax_alpha[i,t,]), col = cols1, clockwise=TRUE,
        main = paste0("True\nCommunity Distribution (i=",i,", t=",t,")"))
    pie(sort(dynbiLCM_Object$softmax_alpha[i,t,]), order(dynbiLCM_Object$softmax_alpha[i,t,]), col = cols2, clockwise=TRUE,
        main = paste0("Estimated\nCommunity Distribution (i=",i,", t=",t,")"))
  } else {
    pie(sort(simulated_data$beta[i,t,]), order(simulated_data$beta[i,t,]), col = cols1, clockwise=TRUE,
        main = paste0("True\nCommunity Distribution (j=",i,", t=",t,")"))
    pie(sort(dynbiLCM_Object$beta[i,t,]), order(dynbiLCM_Object$beta[i,t,]), col = cols2, clockwise=TRUE,
        main = paste0("Estimated\nCommunity Distribution (j=",i,", t=",t,")"))
  }
  
  par(mfrow=c(1,1))
  
}


# # Example -----------------------------------------------------------------
# 
# true_param = random.dynamic.biLCM()
# res = dynamic.biLCM(A = true_param$A)
# plot.compare.dynbiLCM(res, true_param, i = 1, t = 1, group1 = T)

