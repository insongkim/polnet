// Latent space network model.
//
// Implements the model
//   A_{ij} ~ Binom(n, mu_{ij})
// where
//   mu_{ij} = logit(a_i + b_j - ||x_i - y_j||_2^2)
//           = exp(aa_i + ba_j + 2<x_i, y_j>)/(1 + exp(aa_i + ba_j + 2<x_i, y_j>))
//
// The names of the variables in the model are:
//   aa_i : row_factor_adj
//   ba_j : col_factor_adj
//   x_i  : row_embedding
//   y_j  : col_embedding


data {
  int<lower=1> D;     // dimension of latent space
  int<lower=0> n; // number of trials

  int<lower=2> N_row; // number of rows
  int<lower=2> N_col; // number of columns
  int N_fixed_row;    // number of pinned rows
  int N_fixed_col;    // number of pinned columns

  int<lower=0> edges[N_row, N_col];  // connection strength data

  int fixed_row_index[N_fixed_row];             // indices of fixed row agents
  matrix[N_fixed_row, D] fixed_row_embedding;   // positions of fixed row agents

  int fixed_col_index[N_fixed_col];             // indices of fixed column agents
  matrix[D, N_fixed_col] fixed_col_embedding;   // positions of fixed column agents
}

transformed data {
  int flat_ix;
  int flat_edges[N_row * N_col];

  flat_ix = 1;
  for (j in 1:N_col) {
    for (i in 1:(N_row)) {
      flat_edges[flat_ix] = edges[i][j];
      flat_ix = flat_ix + 1;
    }
  }
}

parameters {
  vector<lower=0.01>[D] cov_embedding_diag;

  matrix[N_row, D] row_embedding;
  matrix[D, N_col] col_embedding;

  real mu_col_factor_adj;
  real<lower=0.01> var_row_factor_adj;
  real<lower=0.01> var_col_factor_adj;

  vector[N_row] row_factor_adj;
  row_vector[N_col] col_factor_adj;
}

model {
  vector[N_row * N_col] means;
  int fixed_row_flag;
  int fixed_col_flag;

  for (i in 1:N_row) {
    fixed_row_flag = 0;
    for (j in 1:N_fixed_row) {
      if (i == fixed_row_index[j]) {
        row_embedding[i,] ~ normal(fixed_row_embedding[j,], 1e-4);
        fixed_row_flag = 1;
      }
    }
    if (fixed_row_flag == 0) {
      row_embedding[i,] ~ normal(0.0, cov_embedding_diag);
    }
  }

  for (j in 1:N_col) {
    fixed_col_flag = 0;
    for (k in 1:N_fixed_col) {
      if (j == fixed_col_index[k]) {
        col_embedding[,j] ~ normal(fixed_col_embedding[,k], 1e-4);
        fixed_col_flag = 1;
      }
    }
    if (fixed_col_flag == 0) {
      col_embedding[,j] ~ normal(0.0, cov_embedding_diag);
    }
  }

  row_factor_adj ~ normal(0.0, var_row_factor_adj);
  col_factor_adj ~ normal(mu_col_factor_adj, var_col_factor_adj);

  means = to_vector(
    rep_matrix(row_factor_adj, N_col) +
    rep_matrix(col_factor_adj, N_row) +
    2.0 * row_embedding * col_embedding);

  flat_edges ~ binomial_logit(means);
}
