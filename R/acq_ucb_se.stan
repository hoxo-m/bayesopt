functions {
  real kernel_gp_squared_exponential(row_vector x1, row_vector x2, vector theta) {
    vector[cols(x1)] x1a;
    vector[cols(x2)] x2a;
    real value;
    x1a = to_vector(x1) ./ theta;
    x2a = to_vector(x2) ./ theta;
    value = exp(- squared_distance(x1a, x2a) );
    return if_else(value < 0, 0.0, value);
  }
  matrix compute_covariance_matrix(matrix x_mat, real nu, real theta0, vector theta) {
    matrix[rows(x_mat), rows(x_mat)] cov_mat;
    real value;
    for (i in 1:rows(x_mat)) {
      for (j in i:rows(x_mat)) {
        if (i == j) {
          value = theta0 + nu;
          cov_mat[i, j] = if_else(value < 0 , 0.0, value);
        } else {
          value = theta0 * kernel_gp_squared_exponential(x_mat[i,], x_mat[j, ], theta);
          cov_mat[i, j] = if_else(value < 0 , 0.0, value);
          cov_mat[j, i] = if_else(value < 0 , 0.0, value);
        }
      }
    }
    return cov_mat;
  }
  vector acq_ucb_se(matrix new_x_mat, matrix x_mat, vector y, real m, real nu, real theta0, vector theta, real kappa) {
    matrix[rows(x_mat), rows(x_mat)] cov_mat;
    matrix[rows(x_mat), rows(x_mat)] inv_cov_mat;
    vector[rows(y)] y_colvec;
    matrix[rows(new_x_mat), rows(x_mat)] k;
    vector[rows(new_x_mat)] mu;
    vector[rows(new_x_mat)] sigma;
    real value;
    cov_mat = compute_covariance_matrix(x_mat, nu, theta0, theta);
    inv_cov_mat = inverse_spd(cov_mat);
    y_colvec = y - m;
    for (i in 1:rows(new_x_mat)) {
      for (j in 1:rows(x_mat)) {
        k[i, j] = kernel_gp_squared_exponential(new_x_mat[i, ], x_mat[j, ], theta);
      }
    }
    mu = k * inv_cov_mat * y_colvec;
    for (i in 1:rows(k)) {
      value = theta0 - quad_form(inv_cov_mat, to_vector(k[i,]));
      sigma[i] = if_else(value < 0, 0, sqrt(value));
    }
    return mu + kappa * sigma;
  }
}
data {
  int<lower=0> n_sample;
  int<lower=0> n_dim;
  int<lower=0> n_sample_new;
  matrix[n_sample, n_dim] x_mat;
  vector[n_sample] y;
  matrix[n_sample_new, n_dim] new_x_mat;
  real<lower=0> kappa;
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
generated quantities {
  vector[rows(new_x_mat)] acq;
  acq = acq_ucb_se(new_x_mat, x_mat, y, mu, nu, theta0, theta, kappa);
}
