data {
int<lower=2> N_row; // number of rows
int<lower=2> N_col; // number of columns
int<lower=1> D; // dimensionality of latent space
int<lower=0> edges[N_row, N_col]; // connection strength data
}
transformed data {
int flat_ix;
int flat_edges[N_row * N_col];
flat_ix = 1;
for (i in 1:N_row) {
for (j in 1:N_col) {
flat_edges[flat_ix] = edges[i][j];
flat_ix = flat_ix + 1;
}
}
}
parameters {
vector[D] mu_col_embedding;
vector<lower=0.05>[D] cov_row_embedding_diag;
cov_matrix[D] cov_col_embedding;
row_vector[D] row_embedding[N_row];
row_vector[D] col_embedding[N_col];
real mu_row_factor;
real<lower=0.05> sigma_row_factor;
real<lower=0.05> sigma_col_factor;
vector[N_row] row_factor;
vector[N_col] col_factor;
}
model {
int flat_jx;
vector[N_row * N_col] flat_log_means;
row_factor ~ normal(mu_row_factor, sigma_row_factor);
col_factor ~ normal(0.0, sigma_col_factor);
row_embedding ~ multi_normal(
rep_vector(0.0, D), diag_matrix(cov_row_embedding_diag));
col_embedding ~ multi_normal(
mu_col_embedding, cov_col_embedding);
flat_jx = 1;
for (i in 1:N_row) {
for (j in 1:N_col) {
flat_log_means[flat_jx] =
row_factor[i] + col_factor[j] - distance(
row_embedding[i], col_embedding[j]);
flat_jx = flat_jx + 1;
}
}
flat_edges ~ poisson_log(flat_log_means);
}

