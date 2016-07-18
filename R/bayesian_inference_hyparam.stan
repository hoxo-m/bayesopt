functions {
  real kernel_gp_squared_exponential(row_vector x1, row_vector x2, vector theta) {
    vector[cols(x1)] x1a;
    vector[cols(x2)] x2a;
    x1a = to_vector(x1) ./ theta;
    x2a = to_vector(x2) ./ theta;
    return exp(- squared_distance(x1a, x2a) );
  }
  matrix compute_covariance_matrix(matrix x_mat, real nu, real theta0, vector theta) {
    matrix[rows(x_mat), rows(x_mat)] cov_mat;
    for (i in 1:rows(x_mat)) {
      for (j in i:rows(x_mat)) {
        real value;
        value = kernel_gp_squared_exponential(x_mat[i,], x_mat[j, ], theta);
        if (i == j) {
          cov_mat[i, j] = value + nu;
        } else {
          cov_mat[i, j] = value;
          cov_mat[j, i] = value;
        }
      }
    }
    return theta0 * cov_mat;
  }
}
data {
  int<lower=0> n_sample;
  int<lower=0> n_dim;
  matrix[n_sample, n_dim] x_mat;
  vector[n_sample] y;
}
parameters {
  real mu;
  real<lower=0> nu;
  real<lower=0> theta0;
  vector<lower=0>[n_dim] theta;
}
model {
  matrix[n_sample, n_sample] cov_mat;
  cov_mat = compute_covariance_matrix(x_mat, nu, theta0, theta);
  y ~ multi_normal(rep_vector(mu, n_sample), cov_mat);
  mu ~ cauchy(0, 5);
  nu ~ cauchy(0, 5);
  theta0 ~ cauchy(0, 5);
  theta ~ cauchy(0, 5);
}
