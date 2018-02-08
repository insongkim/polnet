#'@param trained_obj A trained object of class LSNM
#'@param client_space A matrix representing the true latent client space. This matrix should have rows equal to the number of clients, and columns equal to the dimensionality of the latent space, D. 
#'@param legislator_space A matrix representing the true latent legislator space. This matrix should have rows equal to the number of legislators, and columns equal to the dimensionality of the latent space, D. 
#'@return Does not return an object. Prints the proportion of latent space estimates that fell within the credible interval as well as the average error from the true latent space estimates. 
#'

#'@useDynLib polnet, .registration = TRUE
#'@export


comparison <- function(trained_obj, client_space, legislator_space){
  lsnmobj = polnet::get_latent_space(trained_obj)
  l_ordered = lsnmobj[order(rownames(lsnmobj)), ]
  tru_pars = c(as.vector(t(client_space)), as.vector(t(legislator_space)))
  perc_in_cred = sum(tru_pars < l_ordered$`90%` & tru_pars > l_ordered$`10%`) / length(tru_pars)
  ave_marg_error = mean(abs(l_ordered$Mean - tru_pars))
  paste0('A proportion of ', perc_in_cred, 'of latent space parameters fell within their credible interval for an average error of ', ave_marg_error)
}