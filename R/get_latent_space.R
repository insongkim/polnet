#'@param LSNM_Object A trained object of class LSNM
#'@param low_perc The bottom range of the desired credible interval, defaults to 0.1
#'@param high_perc The top range of the credible interval, defaults to 0.9
#'@return A matrix that includes the mean, standard deviation, and credible interval of the latent space estimated by the LSNM algorithm. The row embeddings are the client latent space positions, while the column embeddings are the legislator latent space positions. 
#'

#'@useDynLib polnet, .registration = TRUE
#'@export

get_latent_space <- function(LSNM_Object, low_perc = 0.1, high_perc = 0.9){
  df_fit = as.data.frame(LSNM_Object$stan_fitted_model)
  nms = df_fit[ , grepl( "^col_embedding|^row_embedding" , names(df_fit) )]
  final_df = as.data.frame(cbind(colMeans(nms), apply(nms, 2, sd), apply(nms, 2, quantile, low_perc), apply(nms, 2, quantile, high_perc)))
  names(final_df) = c('Mean', 'SD', '10%', '90%')
  return(final_df)
}