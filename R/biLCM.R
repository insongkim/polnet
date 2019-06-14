#'@name biLSM
#'@param A Matrix of connection strength as counts
#'@param m Number of clients
#'@param n Number of politicians
#'@param k Number of link communities
#'@return A list specifying the membership weights, client membership distributions, and politician membership distributions

#'@export

biLCM <- function(A, m, n, k){

loglike_old <- 1
loglike_current <- 0
converge <- .000001

# Initializes variables

q <- array(1, c(m,n,k)) 
kappa <- array(runif(k, min = 1, max = 10), k)
alpha <-  array(runif(m*k, min = 1, max = 10), c(m,k))
for (z in 1:k){
	alpha[, z] <- alpha[, z]/sum(alpha[, z])
}

beta <-  array(runif(n*k, min = 1, max = 10), c(n,k))
for (z in 1:k){
	beta[, z] <- beta[, z]/sum(beta[, z])
}

# Begins iterative procedure, stops when desired convergence has been reached

while (abs(loglike_current - loglike_old) > converge) {
	
	loglike_old <- loglike_current	
	
	# Update q
	for (i in 1:m){
		for (j in 1:n){
			for (z in 1:k){
				q[i, j, z] <- kappa[z]*alpha[i,z]*beta[j,z]
			}
			q[i, j, ] <- q[i, j, ]/sum(q[i, j, ])
		}
		
	}
	
	# Update kappa
	for (z in 1:k){
		numer <- 0
		denom <- 0
		for (i in 1:m) {			
			for (j in 1:n) {
				numer <- numer + A[i, j]*q[i, j, z]
				denom <- denom + alpha[i, z]*beta[j, z]
			}
		} 
		kappa[z] <- numer/denom
	}
	# Update alpha
	for (i in 1:m){
		for (z in 1:k){
			alpha[i, z] <- 0
			for (j in 1:n){
				alpha[i, z] <- alpha[i, z] + A[i, j]*q[i, j, z]
			}
			
		}
	}
	for (z in 1:k){
		alpha[, z] <- alpha[, z]/sum(alpha[, z])
	}
	# Update beta
	for (j in 1:n){
		for (z in 1:k){
			beta[j, z] <- 0
			for (i in 1:m){
				beta[j, z] <- beta[j, z] + A[i, j]*q[i, j, z]
			}
			
		}
	}
	for (z in 1:k){
		beta[, z] <- beta[, z]/sum(beta[, z])
	}
	
	# Calculate log likelihood
	t1 <- 0
	for (i in 1:m) {
		for (j in 1:n) {
			t2 <- 0
			for (z in 1:k) {
				t2 <- t2 + kappa[z]*alpha[i, z]*beta[j, z]
			}
			t1 <- t1+ A[i, j]*t2
		}
	}
	t3 <- 0
	for (i in 1:m) {
		for (j in 1:n) {
			for (z in 1:k) {
				t3 <- t3 + kappa[z]*alpha[i, z]*beta[j, z]
			}
		}
	}
	loglike_current <- t1 - t3
}

	answer <- list(kappa = kappa, alpha = alpha, beta = beta)
	return(answer)
}
