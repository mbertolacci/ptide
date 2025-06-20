functions {
  matrix cov_periodic(
    real[] x,
    real tau,
    real ell,
    real period
  ) {
    int n = size(x);
    matrix[n, n] output;
    for (i in 1 : n) {
      for (j in (i + 1) : n) {
        output[i, j] = square(tau) * exp(
          -2 * square(
            sin(pi() * (x[i] - x[j]) / period)
          ) / square(ell)
        );
        output[j, i] = output[i, j];
      }
      output[i, i] = square(tau);
    }
    return output;
  }

  matrix cov_quasi_periodic(
    real[] x,
    real tau,
    real ell,
    real ell_decay,
    real period
  ) {
    int n = size(x);
    matrix[n, n] output = gp_exp_quad_cov(x, tau, ell_decay);
    for (i in 1 : n) {
      for (j in (i + 1) : n) {
        output[i, j] = output[i, j] * exp(
          -2 * square(
            sin(pi() * (x[i] - x[j]) / period)
          ) / square(ell)
        );
        output[j, i] = output[i, j];
      }
    }
    return output;
  }

  // Solves (LL')X = b for X, via L^{-T} L^{-1} b
  vector chol_solve_L_b(
    matrix L,
    vector b
  ) {
    return mdivide_right_tri_low(mdivide_left_tri_low(L, b)', L)';
  }

  matrix kronecker_product_sym(
    matrix A,
    matrix B
  ) {
    int n_a = rows(A);
    int n_b = rows(B);
    matrix[n_a * n_b, n_a * n_b] output;
    for (col_A in 1 : n_a) {
      for (row_A in col_A : n_a) {
        matrix[n_b, n_b] output_block = A[row_A, col_A] * B;

        output[
          ((col_A - 1) * n_b + 1) : (col_A * n_b),
          ((row_A - 1) * n_b + 1) : (row_A * n_b)
        ] = output_block;
        output[
          ((row_A - 1) * n_b + 1) : (row_A * n_b),
          ((col_A - 1) * n_b + 1) : (col_A * n_b)
        ] = output_block;
      }
    }
    return output;
  }

  matrix compose_covariance(
    vector sigma,
    matrix L
  ) {
    return quad_form_diag(L * L', sigma);
  }

  int get_block_start(int i, int[] block_last_index) {
    if (i == 1) {
      return 1;
    } else {
      return block_last_index[i - 1] + 1;
    }
  }

  int get_n_block_max(int[] block_last_index) {
    int N_blocks = size(block_last_index);
    int current_block_start_d = 1;
    int n_block_max = 0;
    for (i in 1:N_blocks) {
      int N_block_i = block_last_index[i] - current_block_start_d + 1;
      if (N_block_i > n_block_max) {
        n_block_max = N_block_i;
      }
      current_block_start_d = block_last_index[i] + 1;
    }
    return n_block_max;
  }

  matrix process_covariance(
    data real[] x,
    data int d,
    vector[] exponential_sigma,
    matrix[] exponential_L,
    real[] exponential_ell,
    vector[] matern32_sigma,
    matrix[] matern32_L,
    real[] matern32_ell,
    vector[] matern52_sigma,
    matrix[] matern52_L,
    real[] matern52_ell,
    vector[] squared_exponential_sigma,
    matrix[] squared_exponential_L,
    real[] squared_exponential_ell,
    vector[] quasi_periodic_sigma,
    matrix[] quasi_periodic_L,
    real[] quasi_periodic_ell_periodic,
    real[] quasi_periodic_ell_decay,
    vector[] white_noise_sigma,
    matrix[] white_noise_L,
    data real[] exponential_length_dilation,
    data real[] matern32_length_dilation,
    data real[] matern52_length_dilation,
    data real[] squared_exponential_length_dilation,
    data real[] quasi_periodic_period,
    data real[] quasi_periodic_length_dilation
  ) {
    int n = size(x);
    int n_stacked = n * d;
    int n_exponential_kernels = size(exponential_ell);
    int n_matern32_kernels = size(matern32_ell);
    int n_matern52_kernels = size(matern52_ell);
    int n_squared_exponential_kernels = size(squared_exponential_ell);
    int n_quasi_periodic_kernels = size(quasi_periodic_ell_periodic);
    int n_white_noise_kernels = size(white_noise_sigma);
    matrix[n_stacked, n_stacked] output = rep_matrix(0, n_stacked, n_stacked);
    matrix[d, d] noise_Sigma;

    if (n_exponential_kernels > 0) {
      for (i in 1 : n_exponential_kernels) {
        output += kronecker_product_sym(
          compose_covariance(
            exponential_sigma[i],
            exponential_L[i]
          ),
          gp_exponential_cov(
            x,
            1,
            exponential_ell[i] / exponential_length_dilation[i]
          )
        );
      }
    }

    if (n_matern32_kernels > 0) {
      for (i in 1 : n_matern32_kernels) {
        output += kronecker_product_sym(
          compose_covariance(
            matern32_sigma[i],
            matern32_L[i]
          ),
          gp_matern32_cov(
            x,
            1,
            matern32_ell[i] / matern32_length_dilation[i]
          )
        );
      }
    }
    if (n_matern52_kernels > 0) {
      for (i in 1 : n_matern52_kernels) {
        output += kronecker_product_sym(
          compose_covariance(
            matern52_sigma[i],
            matern52_L[i]
          ),
          gp_matern52_cov(
            x,
            1,
            matern52_ell[i] / matern52_length_dilation[i]
          )
        );
      }
    }

    if (n_squared_exponential_kernels > 0) {
      for (i in 1 : n_squared_exponential_kernels) {
        output += kronecker_product_sym(
          compose_covariance(
            squared_exponential_sigma[i],
            squared_exponential_L[i]
          ),
          gp_exp_quad_cov(
            x,
            1,
            squared_exponential_ell[i] / squared_exponential_length_dilation[i]
          )
        );
      }
    }

    if (n_quasi_periodic_kernels > 0) {
      for (i in 1 : n_quasi_periodic_kernels) {
        output += kronecker_product_sym(
          compose_covariance(quasi_periodic_sigma[i], quasi_periodic_L[i]),
          cov_quasi_periodic(
            x,
            1,
            quasi_periodic_ell_periodic[i],
            quasi_periodic_ell_decay[i] / quasi_periodic_length_dilation[i],
            quasi_periodic_period[i]
          )
        );
      }
    }

    if (n_white_noise_kernels > 0) {
      for (i in 1 : n_white_noise_kernels) {
        noise_Sigma = compose_covariance(white_noise_sigma[i], white_noise_L[i]);
        for (row_d in 1 : d) {
          for (col_d in 1 : d) {
            for (j in 1 : n) {
              output[
                (row_d - 1) * n + j,
                (col_d - 1) * n + j
              ] += noise_Sigma[row_d, col_d];
            }
          }
        }
      }
    }

    return output;
  }

  vector vecchia_reduce_sum(
    int start,
    int end,
    int grain_size,
    int strategy,
    vector[] exponential_sigma,
    matrix[] exponential_L,
    real[] exponential_ell,
    vector[] matern32_sigma,
    matrix[] matern32_L,
    real[] matern32_ell,
    vector[] matern52_sigma,
    matrix[] matern52_L,
    real[] matern52_ell,
    vector[] squared_exponential_sigma,
    matrix[] squared_exponential_L,
    real[] squared_exponential_ell,
    vector[] quasi_periodic_sigma,
    matrix[] quasi_periodic_L,
    real[] quasi_periodic_ell_periodic,
    real[] quasi_periodic_ell_decay,
    vector[] white_noise_sigma,
    matrix[] white_noise_L,
    data real[] x,
    data int d,
    data vector y_stacked,
    data vector y_stacked_tilde,
    data matrix X_stacked,
    data real[] exponential_length_dilation,
    data real[] matern32_length_dilation,
    data real[] matern52_length_dilation,
    data real[] squared_exponential_length_dilation,
    data real[] quasi_periodic_period,
    data real[] quasi_periodic_length_dilation,
    data int[] block_indices,
    data int[] block_last_index,
    data int[] block_n_responses,
    data int[] block_kernel_group
  );

  vector vecchia_partial_sums(
    int start,
    int end,
    vector[] exponential_sigma,
    matrix[] exponential_L,
    real[] exponential_ell,
    vector[] matern32_sigma,
    matrix[] matern32_L,
    real[] matern32_ell,
    vector[] matern52_sigma,
    matrix[] matern52_L,
    real[] matern52_ell,
    vector[] squared_exponential_sigma,
    matrix[] squared_exponential_L,
    real[] squared_exponential_ell,
    vector[] quasi_periodic_sigma,
    matrix[] quasi_periodic_L,
    real[] quasi_periodic_ell_periodic,
    real[] quasi_periodic_ell_decay,
    vector[] white_noise_sigma,
    matrix[] white_noise_L,
    data real[] x,
    data int d,
    data vector y_stacked,
    data vector y_stacked_tilde,
    data matrix X_stacked,
    data real[] exponential_length_dilation,
    data real[] matern32_length_dilation,
    data real[] matern52_length_dilation,
    data real[] squared_exponential_length_dilation,
    data real[] quasi_periodic_period,
    data real[] quasi_periodic_length_dilation,
    data int[] block_indices,
    data int[] block_last_index,
    data int[] block_n_responses,
    data int[] block_kernel_group
  ) {
    int n = size(x);
    int p_stacked = cols(X_stacked);
    int n_kernel_groups = max(block_kernel_group[start:end]);
    int group_n_blocks[n_kernel_groups] = rep_array(0, n_kernel_groups);
    int group_block_size[n_kernel_groups];
    int group_blocks_processed[n_kernel_groups] = rep_array(0, n_kernel_groups);
    int group_block_indices[end - start + 1];
    // Accumulated outputs
    real log_det = 0;
    real y_tildet_Q_y_tilde = 0;
    matrix[p_stacked, p_stacked] Xt_Q_X = rep_matrix(0, p_stacked, p_stacked);
    vector[p_stacked] Xt_Q_y = rep_vector(0, p_stacked);
    vector[p_stacked] Xt_Q_y_tilde = rep_vector(0, p_stacked);

    for (j in start : end) {
      group_n_blocks[block_kernel_group[j]] += 1;
      group_block_size[block_kernel_group[j]] = block_last_index[j] - get_block_start(j, block_last_index) + 1;
    }

    for (j in start : end) {
      int group = block_kernel_group[j];
      int group_start = 1;
      if (group > 1) {
        group_start = sum(group_n_blocks[1 : (group - 1)]) + 1;
      }
      group_block_indices[
        group_start + group_blocks_processed[group]
      ] = j;
      group_blocks_processed[group] += 1;
    }

    for (i in 1 : n_kernel_groups) {
      if (group_n_blocks[i] > 0) {
        int current_block_size = group_block_size[i];
        int current_block_size_stacked = d * current_block_size;
        int current_n_blocks = group_n_blocks[i];

        matrix[current_block_size_stacked, current_block_size_stacked] K_group;
        matrix[current_block_size_stacked, current_block_size_stacked] L_group;

        int current_block_indices[current_block_size];
        // First p_stacked * current_n_blocks columns come from X, second
        // current_n_blocks columns come from y, and last current_n_blocks
        // come from y_tilde
        matrix[
          current_block_size_stacked,
          p_stacked * current_n_blocks + 2 * current_n_blocks
        ] Z_group;
        matrix[
          current_block_size_stacked,
          p_stacked * current_n_blocks + 2 * current_n_blocks
        ] L_inv_Z_group;

        for (j in 1 : group_n_blocks[i]) {
          int block_index = group_block_indices[
            sum(group_n_blocks[1 : (i - 1)]) + j
          ];
          int block_start = get_block_start(block_index, block_last_index);
          int current_block_indices_stacked[current_block_size_stacked];
          current_block_indices = block_indices[
            block_start:block_last_index[block_index]
          ];

          for (k in 1 : d) {
            for (l in 1 : current_block_size) {
              current_block_indices_stacked[(k - 1) * current_block_size + l] = (
                (k - 1) * n + current_block_indices[l]
              );
            }
          }

          Z_group[
            :,
            ((j - 1) * p_stacked + 1) : (j * p_stacked)
          ] = X_stacked[current_block_indices_stacked, :];
          Z_group[
            :,
            current_n_blocks * p_stacked + j
          ] = y_stacked[current_block_indices_stacked];
          Z_group[
            :,
            current_n_blocks * p_stacked + current_n_blocks + j
          ] = y_stacked_tilde[current_block_indices_stacked];
        }

        K_group = process_covariance(
          // NOTE: uses the last entry to determine x, since it's all the same
          x[current_block_indices],
          d,
          exponential_sigma,
          exponential_L,
          exponential_ell,
          matern32_sigma,
          matern32_L,
          matern32_ell,
          matern52_sigma,
          matern52_L,
          matern52_ell,
          squared_exponential_sigma,
          squared_exponential_L,
          squared_exponential_ell,
          quasi_periodic_sigma,
          quasi_periodic_L,
          quasi_periodic_ell_periodic,
          quasi_periodic_ell_decay,
          white_noise_sigma,
          white_noise_L,
          exponential_length_dilation,
          matern32_length_dilation,
          matern52_length_dilation,
          squared_exponential_length_dilation,
          quasi_periodic_period,
          quasi_periodic_length_dilation
        );
        L_group = cholesky_decompose(K_group);
        L_inv_Z_group = mdivide_left_tri_low(L_group, Z_group);

        for (j in 1 : group_n_blocks[i]) {
          int block_index = group_block_indices[
            sum(group_n_blocks[1 : (i - 1)]) + j
          ];
          int block_start = get_block_start(block_index, block_last_index);
          int n_parents_current_block = current_block_size - block_n_responses[block_index];
          int response_rows_current_block_stacked[d * block_n_responses[block_index]];
          current_block_indices = block_indices[
            block_start:block_last_index[block_index]
          ];

          for (k in 1 : d) {
            for (l in 1 : block_n_responses[block_index]) {
              response_rows_current_block_stacked[
                (k - 1) * block_n_responses[block_index] + l
              ] = (k - 1) * current_block_size + n_parents_current_block + l;
            }
          }

          log_det += (
            2 * sum(log(diagonal(L_group[
              response_rows_current_block_stacked,
              response_rows_current_block_stacked
            ])))
          );
          y_tildet_Q_y_tilde += (
            sum(square(L_inv_Z_group[
              response_rows_current_block_stacked,
              current_n_blocks * p_stacked + current_n_blocks + j
            ]))
          );
          Xt_Q_y += (
            L_inv_Z_group[
              response_rows_current_block_stacked,
              ((j - 1) * p_stacked + 1) : (j * p_stacked)
            ]' * L_inv_Z_group[
              response_rows_current_block_stacked,
              current_n_blocks * p_stacked + j
            ]
          );
          Xt_Q_y_tilde += (
            L_inv_Z_group[
              response_rows_current_block_stacked,
              ((j - 1) * p_stacked + 1) : (j * p_stacked)
            ]' * L_inv_Z_group[
              response_rows_current_block_stacked,
              current_n_blocks * p_stacked + current_n_blocks + j
            ]
          );
          Xt_Q_X += crossprod(L_inv_Z_group[
            response_rows_current_block_stacked,
            ((j - 1) * p_stacked + 1) : (j * p_stacked)
          ]);
        }
      }
    }

    vector[2 + p_stacked * p_stacked + 2 * p_stacked] output;
    output[1] = log_det;
    output[2] = y_tildet_Q_y_tilde;
    output[3:(2 + p_stacked * p_stacked)] = to_vector(Xt_Q_X);
    output[(3 + p_stacked * p_stacked):(2 + p_stacked * p_stacked + p_stacked)] = to_vector(Xt_Q_y);
    output[(3 + p_stacked * p_stacked + p_stacked):] = to_vector(Xt_Q_y_tilde);
    return output;
  }
}
data {
  int<lower=1> n;
  int<lower=1> p;
  int<lower=1> d;
  real x[n];
  matrix[n, d] y;
  matrix[n, p] X;

  vector[p * d] beta_prior_mean;
  matrix[p * d, p * d] beta_prior_precision;

  int<lower=0> n_exponential_kernels;
  real<lower=0> exponential_length_dilation[n_exponential_kernels];
  real<lower=0> exponential_sigma_squared_a[n_exponential_kernels];
  real<lower=0> exponential_sigma_squared_b[n_exponential_kernels];
  real<lower=0> exponential_L_shape[n_exponential_kernels];
  real<lower=0> exponential_ell_a[n_exponential_kernels];
  real<lower=0> exponential_ell_b[n_exponential_kernels];

  int<lower=0> n_matern32_kernels;
  real<lower=0> matern32_length_dilation[n_matern32_kernels];
  real<lower=0> matern32_sigma_squared_a[n_matern32_kernels];
  real<lower=0> matern32_sigma_squared_b[n_matern32_kernels];
  real<lower=0> matern32_L_shape[n_matern32_kernels];
  real<lower=0> matern32_ell_a[n_matern32_kernels];
  real<lower=0> matern32_ell_b[n_matern32_kernels];

  int<lower=0> n_matern52_kernels;
  real<lower=0> matern52_length_dilation[n_matern52_kernels];
  real<lower=0> matern52_sigma_squared_a[n_matern52_kernels];
  real<lower=0> matern52_sigma_squared_b[n_matern52_kernels];
  real<lower=0> matern52_L_shape[n_matern52_kernels];
  real<lower=0> matern52_ell_a[n_matern52_kernels];
  real<lower=0> matern52_ell_b[n_matern52_kernels];

  int<lower=0> n_squared_exponential_kernels;
  real<lower=0> squared_exponential_length_dilation[n_squared_exponential_kernels];
  real<lower=0> squared_exponential_sigma_squared_a[n_squared_exponential_kernels];
  real<lower=0> squared_exponential_sigma_squared_b[n_squared_exponential_kernels];
  real<lower=0> squared_exponential_L_shape[n_squared_exponential_kernels];
  real<lower=0> squared_exponential_ell_a[n_squared_exponential_kernels];
  real<lower=0> squared_exponential_ell_b[n_squared_exponential_kernels];

  int<lower=0> n_quasi_periodic_kernels;
  real<lower=0> quasi_periodic_period[n_quasi_periodic_kernels];
  real<lower=0> quasi_periodic_length_dilation[n_quasi_periodic_kernels];
  real<lower=0> quasi_periodic_sigma_squared_a[n_quasi_periodic_kernels];
  real<lower=0> quasi_periodic_sigma_squared_b[n_quasi_periodic_kernels];
  real<lower=0> quasi_periodic_L_shape[n_quasi_periodic_kernels];
  real<lower=0> quasi_periodic_ell_periodic_a[n_quasi_periodic_kernels];
  real<lower=0> quasi_periodic_ell_periodic_b[n_quasi_periodic_kernels];
  real<lower=0> quasi_periodic_ell_decay_a[n_quasi_periodic_kernels];
  real<lower=0> quasi_periodic_ell_decay_b[n_quasi_periodic_kernels];

  int<lower=0> n_white_noise_kernels;
  real<lower=0> white_noise_sigma_squared_a[n_white_noise_kernels];
  real<lower=0> white_noise_sigma_squared_b[n_white_noise_kernels];
  real<lower=0> white_noise_L_shape[n_white_noise_kernels];

  // Vecchia approximation blocks
  int<lower=1> n_indices;
  int<lower=1> n_blocks;
  int block_indices[n_indices];
  int block_last_index[n_blocks];
  int block_n_responses[n_blocks];
  int block_kernel_group[n_blocks];

  // Computation strategy
  int use_parallel;
  int<lower=1> grain_size;
  int strategy;
}
transformed data {
  int n_block_max = get_n_block_max(block_last_index);
  int n_stacked = n * d;
  int p_stacked = p * d;
  vector[n_stacked] y_stacked = to_vector(y);
  vector[n_stacked] y_stacked_tilde;
  matrix[n_stacked, p_stacked] X_stacked;
  matrix[p_stacked, p_stacked] L_beta_prior_precision = cholesky_decompose(
    beta_prior_precision
  );

  X_stacked = rep_matrix(0, n_stacked, p_stacked);
  for (i in 1 : d) {
    for (j in 1 : p) {
      int stacked_row_start = (i - 1) * n + 1;
      int stacked_col = (i - 1) * p + j;
      X_stacked[
        stacked_row_start : (stacked_row_start + n - 1),
        stacked_col
      ] = X[ , j];
    }
  }

  y_stacked_tilde = y_stacked - X_stacked * beta_prior_mean;
}
parameters {
  vector<lower=0>[d] exponential_sigma_squared[n_exponential_kernels];
  cholesky_factor_corr[d] exponential_L[n_exponential_kernels];
  real<lower=0> exponential_ell[n_exponential_kernels];

  vector<lower=0>[d] matern32_sigma_squared[n_matern32_kernels];
  cholesky_factor_corr[d] matern32_L[n_matern32_kernels];
  real<lower=0> matern32_ell[n_matern32_kernels];

  vector<lower=0>[d] matern52_sigma_squared[n_matern52_kernels];
  cholesky_factor_corr[d] matern52_L[n_matern52_kernels];
  real<lower=0> matern52_ell[n_matern52_kernels];

  vector<lower=0>[d] squared_exponential_sigma_squared[n_squared_exponential_kernels];
  cholesky_factor_corr[d] squared_exponential_L[n_squared_exponential_kernels];
  real<lower=0> squared_exponential_ell[n_squared_exponential_kernels];

  vector<lower=0>[d] quasi_periodic_sigma_squared[n_quasi_periodic_kernels];
  cholesky_factor_corr[d] quasi_periodic_L[n_quasi_periodic_kernels];
  real<lower=0> quasi_periodic_ell_periodic[n_quasi_periodic_kernels];
  real<lower=0> quasi_periodic_ell_decay[n_quasi_periodic_kernels];

  vector<lower=0>[d] white_noise_sigma_squared[n_white_noise_kernels];
  cholesky_factor_corr[d] white_noise_L[n_white_noise_kernels];
}
transformed parameters {
  vector[p_stacked] beta_hat;
  real log_marginal;
  {
    vector[d] exponential_sigma[n_exponential_kernels] = sqrt(exponential_sigma_squared);
    vector[d] matern32_sigma[n_matern32_kernels] = sqrt(matern32_sigma_squared);
    vector[d] matern52_sigma[n_matern52_kernels] = sqrt(matern52_sigma_squared);
    vector[d] squared_exponential_sigma[n_squared_exponential_kernels] = sqrt(squared_exponential_sigma_squared);
    vector[d] quasi_periodic_sigma[n_quasi_periodic_kernels] = sqrt(quasi_periodic_sigma_squared);
    vector[d] white_noise_sigma[n_white_noise_kernels] = sqrt(white_noise_sigma_squared);

    vector[
      2 + p_stacked * p_stacked + 2 * p_stacked
    ] components;
    real log_det = 0;
    real y_tildet_Q_y_tilde = 0;
    matrix[p_stacked, p_stacked] Xt_Q_X = rep_matrix(0, p_stacked, p_stacked);
    vector[p_stacked] Xt_Q_y = rep_vector(0, p_stacked);
    vector[p_stacked] Xt_Q_y_tilde = rep_vector(0, p_stacked);
    matrix[p_stacked, p_stacked] L_beta_precision;
    matrix[p_stacked, p_stacked] beta_precision;
    vector[p_stacked] beta_rhs;

    if (use_parallel > 0) {
      components = vecchia_reduce_sum(
        1,
        n_blocks,
        grain_size,
        strategy,
        exponential_sigma,
        exponential_L,
        exponential_ell,
        matern32_sigma,
        matern32_L,
        matern32_ell,
        matern52_sigma,
        matern52_L,
        matern52_ell,
        squared_exponential_sigma,
        squared_exponential_L,
        squared_exponential_ell,
        quasi_periodic_sigma,
        quasi_periodic_L,
        quasi_periodic_ell_periodic,
        quasi_periodic_ell_decay,
        white_noise_sigma,
        white_noise_L,
        x,
        d,
        y_stacked,
        y_stacked_tilde,
        X_stacked,
        exponential_length_dilation,
        matern32_length_dilation,
        matern52_length_dilation,
        squared_exponential_length_dilation,
        quasi_periodic_period,
        quasi_periodic_length_dilation,
        block_indices,
        block_last_index,
        block_n_responses,
        block_kernel_group
      );
    } else {
      components = vecchia_partial_sums(
        1,
        n_blocks,
        exponential_sigma,
        exponential_L,
        exponential_ell,
        matern32_sigma,
        matern32_L,
        matern32_ell,
        matern52_sigma,
        matern52_L,
        matern52_ell,
        squared_exponential_sigma,
        squared_exponential_L,
        squared_exponential_ell,
        quasi_periodic_sigma,
        quasi_periodic_L,
        quasi_periodic_ell_periodic,
        quasi_periodic_ell_decay,
        white_noise_sigma,
        white_noise_L,
        x,
        d,
        y_stacked,
        y_stacked_tilde,
        X_stacked,
        exponential_length_dilation,
        matern32_length_dilation,
        matern52_length_dilation,
        squared_exponential_length_dilation,
        quasi_periodic_period,
        quasi_periodic_length_dilation,
        block_indices,
        block_last_index,
        block_n_responses,
        block_kernel_group
      );
    }

    log_det = components[1];
    y_tildet_Q_y_tilde = components[2];
    Xt_Q_X = to_matrix(
      components[3:(2 + p_stacked * p_stacked)],
      p_stacked,
      p_stacked
    );
    Xt_Q_y = components[
      (3 + p_stacked * p_stacked):(
        2 + p_stacked * p_stacked + p_stacked
      )
    ];
    Xt_Q_y_tilde = components[
      (3 + p_stacked * p_stacked + p_stacked):
    ];
    beta_precision = Xt_Q_X + beta_prior_precision;
    L_beta_precision = cholesky_decompose(beta_precision);
    beta_rhs = Xt_Q_y + beta_prior_precision * beta_prior_mean;
    beta_hat = chol_solve_L_b(L_beta_precision, beta_rhs);
    // This uses the Woodbury matrix identity
    log_marginal = -0.5 * (
      2 * sum(log(diagonal(L_beta_precision)))
      - 2 * sum(log(diagonal(L_beta_prior_precision)))
      + log_det
      + y_tildet_Q_y_tilde
      - sum(square(mdivide_left_tri_low(
        L_beta_precision,
        Xt_Q_y_tilde
      )))
    );
  }
}
model {
  target += log_marginal;

  if (n_exponential_kernels > 0) {
    for (i in 1 : n_exponential_kernels) {
      exponential_sigma_squared[i] ~ inv_gamma(
        exponential_sigma_squared_a[i],
        exponential_sigma_squared_b[i]
      );
      exponential_L[i] ~ lkj_corr_cholesky(
        exponential_L_shape[i]
      );
      exponential_ell[i] ~ inv_gamma(
        exponential_ell_a[i],
        exponential_ell_b[i]
      );
    }
  }

  if (n_matern32_kernels > 0) {
    for (i in 1 : n_matern32_kernels) {
      matern32_sigma_squared[i] ~ inv_gamma(
        matern32_sigma_squared_a[i],
        matern32_sigma_squared_b[i]
      );
      matern32_L[i] ~ lkj_corr_cholesky(
        matern32_L_shape[i]
      );
      matern32_ell[i] ~ inv_gamma(
        matern32_ell_a[i],
        matern32_ell_b[i]
      );
    }
  }

  if (n_matern52_kernels > 0) {
    for (i in 1 : n_matern52_kernels) {
      matern52_sigma_squared[i] ~ inv_gamma(
        matern52_sigma_squared_a[i],
        matern52_sigma_squared_b[i]
      );
      matern52_L[i] ~ lkj_corr_cholesky(
        matern52_L_shape[i]
      );
      matern52_ell[i] ~ inv_gamma(
        matern52_ell_a[i],
        matern52_ell_b[i]
      );
    }
  }

  if (n_squared_exponential_kernels > 0) {
    for (i in 1 : n_squared_exponential_kernels) {
      squared_exponential_sigma_squared[i] ~ inv_gamma(
        squared_exponential_sigma_squared_a[i],
        squared_exponential_sigma_squared_b[i]
      );
      squared_exponential_L[i] ~ lkj_corr_cholesky(
        squared_exponential_L_shape[i]
      );
      squared_exponential_ell[i] ~ inv_gamma(
        squared_exponential_ell_a[i],
        squared_exponential_ell_b[i]
      );
    }
  }

  if (n_quasi_periodic_kernels > 0) {
    for (i in 1 : n_quasi_periodic_kernels) {
      quasi_periodic_sigma_squared[i] ~ inv_gamma(
        quasi_periodic_sigma_squared_a[i],
        quasi_periodic_sigma_squared_b[i]
      );
      quasi_periodic_L[i] ~ lkj_corr_cholesky(
        quasi_periodic_L_shape[i]
      );
      quasi_periodic_ell_periodic[i] ~ inv_gamma(
        quasi_periodic_ell_periodic_a[i],
        quasi_periodic_ell_periodic_b[i]
      );
      quasi_periodic_ell_decay[i] ~ inv_gamma(
        quasi_periodic_ell_decay_a[i],
        quasi_periodic_ell_decay_b[i]
      );
    }
  }

  if (n_white_noise_kernels > 0) {
    for (i in 1 : n_white_noise_kernels) {
      white_noise_sigma_squared[i] ~ inv_gamma(
        white_noise_sigma_squared_a[i],
        white_noise_sigma_squared_b[i]
      );
      white_noise_L[i] ~ lkj_corr_cholesky(
        white_noise_L_shape[i]
      );
    }
  }
}
