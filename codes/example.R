setwd("~/Google Drive/Sooahn/Network/dynamic/codes")
source("dynamicbiLCM.R")

true_param = random.dynamic.biLCM(# number of community = length of z (+int)
                                  k = 3,
                                  # number of member state = length of i (+int)
                                  m = 20,
                                  # number of time period = length of t (+int; Assmp: \underbar{T_i}=1 and \upperbar{T_i}=s \forall i)
                                  s = 5,
                                  # number of bills for each time period = length of j's (+int; Assmp: n_1 = ... = n_s)
                                  n = 30)
res = dynamic.biLCM(A = true_param$A, k = 3, r= 0.0001)

# Initialization with true parameters
res = dynamic.biLCM(A = true_param$A, k = 3, r= 0.0001, 
                    kappa = true_param$kappa,
                    alpha = true_param$alpha,
                    beta = true_param$beta)
res$loglikelihood
# plot.compare.dynbiLCM(res, true_param, i = 1, t = 1, group1 = T)
A = true_param$A
m = dim(A)[1]
s = dim(A)[3]
n = dim(A)[2]
for (tt in 1:s) {
  for (ii in 1:m) {
    pdf(file = paste0("~/Google Drive/Sooahn/Network/dynamic/figs/synthetic/t",tt,"_i",ii,".pdf"))
    plot.compare.dynbiLCM(res, true_param, i = ii, t = tt, group1 = T)
    dev.off()
  }
}
for (tt in 1:s) {
  for (jj in 1:n) {
    pdf(file = paste0("~/Google Drive/Sooahn/Network/dynamic/figs/synthetic/t",tt,"_j",jj,".pdf"))
    plot.compare.dynbiLCM(res, true_param, i = jj, t = tt, group1 = F)
    dev.off()
  }
}
