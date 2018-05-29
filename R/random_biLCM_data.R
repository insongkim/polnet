#'@param m The number of clients
#'@param n The number of politicians
#'@param k The number of link communities
#'@return A list of the adjacency matrix, the membership vectors for the clients, membership vectors for the legislators, and the base weight kappa

#'@export

random_biLCM_data <- function(m, n, k){

generate <- T

while (generate == T){
A <- array(0, c(m,n))
kappa <- sample(5:30, k, replace = TRUE)

#generating alpha
alpha <-  array(0, c(m,k))
num_m <- rpois(m, k/2) #choose number of communities each client belongs to, k/2 communities on average
for (i in 1:m){
		index_i <- sample(1:k, max(1, min(k, num_m[i])), replace = FALSE) 
		for (index in index_i){
			alpha[i, index] <- rpois(1, 10)
		}
}

for (z in 1:k){
	alpha[, z] <- alpha[, z]/sum(alpha[, z])
}

#generating beta
beta <-  array(0, c(n,k))
num_n <- rpois(n, k/2) #choose number of communities each politician belongs to, k/2 communities on average
for (j in 1:n){
		index_j <- sample(1:k, max(1, min(k, num_n[j])), replace = FALSE) 
		for (index in index_j){
			beta[j, index] <- rpois(1, 10)
		}
}

for (z in 1:k){
	beta[, z] <- beta[, z]/sum(beta[, z])
}

#generate A
for (i in 1:m){
	for (j in 1:n){
		mu <- 0
		for (z in 1:k){
			mu <- mu + kappa[z]*alpha[i, z]*beta[j, z]
		}
		A[i, j] <- rpois(1, mu)
	}
}

generate <- F
#check for rows/columns of all 0
if ((0 %in% colSums(A)) | (0 %in% rowSums(A))){
	generate <- T
}

} #generate loop close

return(list(A, alpha, beta, kappa))
}
