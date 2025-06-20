#' Predict using PTide
#'
#' Perform predictions using a PTide fit. The predictions are made at chosen
#' times and will be conditioned on the observed data and the covariance
#' parameters estimated by \code{\link{ptide_fit}}. The covariate coefficients
#' and the coefficients of the tidal harmonics are re-estimated using the
#' observed data.
#'
#' The computational time required for the prediction depends on the number of
#' observations, the number of predicted time points, the choice of Vecchia
#' approximation (see the \code{vecchia_configuration} parameter below), and
#' whether the observed times are regular. The last point is particularly
#' important: if the observed times are regular, many matrix operations can be
#' precomputed and the prediction will be much faster.
#'
#' @param fit A PTide fit object made using \code{\link{ptide_fit}}.
#' @param predicted_df A data frame containing the times to be predicted, and
#' the values of any covariates used in the fit.
#' @param observed_df A data frame containing the observed data that will be
#' used for the prediction.
#' @param include_sd Whether to include the standard deviation of the
#' prediction in the output.
#' @param include_marginal_covariance Whether to include the marginal
#' covariance of the prediction in the output.
#' @param include_covariance Whether to include the full covariance of the
#' prediction in the output.
#' @param include_white_noise Whether to include the white noise kernel in the
#' prediction variance and the samples. This is only relevant if the white noise
#' kernel was included in the fit.
#' @param n_samples The number of samples to draw from the prediction
#' distribution. If set to 0, no samples will be drawn.
#' @param method The method to use for the prediction. See the methods section
#' below.
#' @param vecchia_configuration A list containing the configuration for the
#' Vecchia approximation. See \code{\link{ptide_fit}} for documentation on this
#' field; in general, larger values are better, but take longer to compute.
#' The \code{n_nearest} entry is used for the standard method when the
#' prediction times are not in the future relative to the observed data.
#' @param cache A cache object satisfying the semantics of the \code{cachem}
#' package, such as \code{\link[cachem]{cache_mem}}. This can be used to save
#' the results of expensive computations between calls to \code{ptide_predict}.
#' If omitted, no cache is used.
#' @param cache_resolution The resolution of the cache, in seconds. This is
#' only relevant if \code{cache} is specified. Covariances over differences
#' smaller than this resolution are considered to be identical.
#'
#' @return A list containing the predicted values. The list will contain the
#' following elements:
#' \itemize{
#'   \item \code{fitted}: The fitted values, which is equal to the mean function
#'   evaluated at the prediction times.
#'   \item \code{mean}: The posterior mean of the prediction distribution.
#'   \item \code{sd}: The posterior standard deviation of the prediction
#'   distribution. This will only be included if \code{include_sd} is set to
#'   \code{TRUE}.
#'   \item \code{marginal_covariance}: The posterior covariance matrix at each
#'   time point that is predicted. This will only be included if
#'   \code{include_marginal_covariance} is set to \code{TRUE}.
#'   \item \code{covariance}: The posterior covariance matrix between all
#'   predicted time points. This will only be included if
#'   \code{include_covariance} is set to \code{TRUE}.
#'   \item \code{samples}: Samples from the posterior distribution. This will
#'   only be included if \code{n_samples} is greater than 0.
#' }
#'
#' @section Methods:
#' Two methods are available:
#' \itemize{
#'   \item \code{standard}: This method corresponds to the standard Vecchia
#'     approximation. This is fast when the prediction times are in the
#'     future, and reasonably fast when some are in the past. It is very
#'     accurate if the number of parents in the Vecchia configuration is high
#'     and when the white noise component of the covariance structure is
#'     relatively small.
#'   \item \code{general}: This method corresponds to the general Vecchia
#'     approximation of Katzfuss et. al (2020). This is much slower than the
#'     standard method but may be more accurate in some situations.
#' }
#'
#' @export
ptide_predict <- function(
  fit,
  predicted_df,
  observed_df,
  include_sd = FALSE,
  include_marginal_covariance = FALSE,
  include_covariance = FALSE,
  include_white_noise = TRUE,
  n_samples = 0,
  method = c('standard', 'general'),
  vecchia_configuration = list(
    n_parents = 256L,
    n_responses = 256L,
    n_nearest = 256L
  ),
  cache,
  cache_resolution = 60L
) {
  method <- match.arg(method)

  logger::log_trace('Converting input data to required format')
  observed_model_data <- .ptide_model_data(fit$model, observed_df)
  predicted_model_data <- .ptide_model_data(
    fit$model,
    predicted_df,
    include_response = FALSE
  )
  observed_model_data_long <- .lengthen_tidal_model_data(
    observed_model_data,
    stack = FALSE
  )
  predicted_model_data_long <- .lengthen_tidal_model_data(
    predicted_model_data,
    stack = FALSE
  )

  is_white_noise_only <- (
    length(fit$covariance_functions) == 1L
    && fit$covariance_functions[[1]]$type == 'white_noise'
  )

  method_call <- if (is_white_noise_only) {
    .ptide_predict_white_noise_only
  } else if (method == 'standard') {
    .ptide_predict_standard
  } else if (method == 'general') {
    .ptide_predict_general
  }

  prediction <- method_call(
    fit,
    observed_model_data,
    predicted_model_data,
    observed_model_data_long,
    predicted_model_data_long,
    include_posterior_covariance = (
      include_sd || include_marginal_covariance || include_covariance
    ),
    n_samples = n_samples,
    vecchia_configuration = vecchia_configuration,
    cache = cache,
    cache_resolution = cache_resolution
  )

  n_predicted <- nrow(predicted_model_data$X)
  d <- observed_model_data$d
  response_names <- colnames(observed_model_data$y)

  needs_covariance <- (
    include_sd || include_marginal_covariance || include_covariance
  )

  if (needs_covariance && include_white_noise) {
    covariance_wn <- matrix(0, nrow = d, ncol = d)
    for (covariance_function in fit$covariance_functions) {
      if (covariance_function$type != 'white_noise') next

      covariance_wn <- covariance_wn + tcrossprod(
        diag(sqrt(covariance_function$sigma_squared))
        %*% covariance_function$L
      )
    }
  }

  if (needs_covariance) {
    covariance <- if (include_white_noise) {
      prediction$covariance + kronecker(diag(n_predicted), covariance_wn)
    } else {
      prediction$covariance
    }

    if (include_marginal_covariance || include_covariance) {
      interleave_to_stack <- (
        rep(d * seq_len(n_predicted), d) + rep((-d + 1L) : 0L, each = n_predicted)
      )
      Sigma_y_predicted_stacked <- covariance[
        interleave_to_stack,
        interleave_to_stack
      ]
    }
  }

  uninterleave_value <- function(input) {
    output <- t(matrix(input, nrow = d))
    colnames(output) <- response_names
    output
  }

  logger::log_trace('Formatting output')
  output <- list(
    fitted = uninterleave_value(prediction$fitted),
    mean = uninterleave_value(prediction$mean)
  )

  if (include_sd) {
    output$sd <- uninterleave_value(sqrt(diag(covariance)))
  }

  if (include_marginal_covariance) {
    output$marginal_covariance <- array(NA, dim = c(n_predicted, d, d))
    for (i in 1 : d) {
      for (j in 1 : d) {
        # Grab the diagonal of the (i, j) block
        output$marginal_covariance[, i, j] <- Sigma_y_predicted_stacked[cbind(
          ((i - 1) * n_predicted + 1) : (i * n_predicted),
          ((j - 1) * n_predicted + 1) : (j * n_predicted)
        )]
      }
    }
    dimnames(output$marginal_covariance)[[2]] <- response_names
    dimnames(output$marginal_covariance)[[3]] <- response_names
  }

  if (include_covariance) {
    output$covariance <- array(
      Sigma_y_predicted_stacked,
      dim = c(n_predicted, d, n_predicted, d)
    )
    dimnames(output$covariance)[[2]] <- response_names
    dimnames(output$covariance)[[4]] <- response_names
  }

  if (n_samples > 0) {
    if (include_white_noise) {
      epsilon <- crossprod(
        chol(covariance_wn),
        matrix(rnorm(n_samples * n_predicted * d), nrow = d)
      )
      prediction$samples <- prediction$samples + matrix(
        epsilon,
        ncol = n_samples
      )
    }

    output$samples <- array(NA, dim = c(n_predicted, d, n_samples))
    for (i in seq_len(n_samples)) {
      output$samples[, , i] <- uninterleave_value(prediction$samples[, i])
    }
  }

  output
}

.ptide_predict_white_noise_only <- function(
  fit,
  observed_model_data,
  predicted_model_data,
  observed_model_data_long,
  predicted_model_data_long,
  include_posterior_covariance,
  n_samples,
  vecchia_configuration,
  cache,
  cache_resolution
) {
  n_predicted <- length(predicted_model_data$x)
  d <- observed_model_data$d

  X_o <- observed_model_data_long$X
  X_p <- predicted_model_data_long$X
  y_o <- observed_model_data_long$y

  L <- (
    diag(sqrt(fit$covariance_functions[[1]]$sigma_squared))
    %*% fit$covariance_functions[[1]]$L
  )
  L_inv_y <- as.vector(forwardsolve(L, matrix(y_o, nrow = d)))
  L_inv_X <- matrix(forwardsolve(L, matrix(X_o, nrow = d)), nrow = nrow(X_o))
  beta_hat <- as.vector(solve(crossprod(L_inv_X),crossprod(L_inv_X, L_inv_y)))

  y_predicted_fitted <- as.vector(X_p %*% beta_hat)
  y_predicted_mean <- y_predicted_fitted

  if (include_posterior_covariance) {
    # NOTE: The covariance is supposed to be exclusive of the white noise
    # component, hence zero
    Sigma_y_predicted <- matrix(
      0,
      nrow = d * n_predicted,
      ncol = d * n_predicted
    )
  }

  if (n_samples > 0) {
    y_samples <- matrix(
      y_predicted_mean,
      nrow = d * n_predicted,
      ncol = n_samples
    )
  }

  list(
    fitted = y_predicted_fitted,
    mean = y_predicted_mean,
    covariance = if (include_posterior_covariance) {
      Sigma_y_predicted
    } else {
      NULL
    },
    samples = if (n_samples > 0) y_samples else NULL
  )
}

.ptide_predict_standard <- function(
  fit,
  observed_model_data,
  predicted_model_data,
  observed_model_data_long,
  predicted_model_data_long,
  include_posterior_covariance,
  n_samples,
  vecchia_configuration,
  cache,
  cache_resolution
) {
  n_parents <- vecchia_configuration$n_parents
  n_responses <- vecchia_configuration$n_responses
  n_nearest <- vecchia_configuration$n_nearest
  d <- observed_model_data$d
  n_observed <- length(observed_model_data$x)
  n_predicted <- length(predicted_model_data$x)
  has_cache <- !missing(cache)

  local_cache <- if (has_cache) {
    cache
  } else {
    cachem::cache_mem(evict = 'fifo')
  }

  # If this is a forecast, we can use a sequential AR-style parent scheme
  # where the forecast times come after the observed times.
  is_forecast <- (max(observed_model_data$x) <= min(predicted_model_data$x))
  stopifnot(!is.unsorted(observed_model_data$x))
  if (is_forecast) {
    stopifnot(!is.unsorted(predicted_model_data$x))
  }

  logger::log_trace('Constructing parents')
  x <- c(observed_model_data$x, predicted_model_data$x)
  if (is_forecast) {
    x_sequential <- x
  } else {
    x_sequential <- observed_model_data$x
  }
  n_sequential_blocks <- ceiling(
    (length(x_sequential) - n_parents) / n_responses
  )
  sequential_block_responses <- c(
    list(1L : (n_parents + n_responses)),
    lapply(2L : n_sequential_blocks, function(i) {
      ((i - 1L) * n_responses + n_parents + 1L) : min(
        i * n_responses + n_parents,
        length(x)
      )
    })
  )
  sequential_block_parents <- c(
    list(integer(0)),
    lapply(2L : n_sequential_blocks, function(i) {
      ((i - 1L) * n_responses + 1L) : ((i - 1L) * n_responses + n_parents)
    })
  )
  if (is_forecast) {
    singleton_responses <- NULL
    singleton_parents <- NULL
  } else {
    predicted_nn <- parent_knn(
      x,
      n_nearest
    )[(n_observed + 1L) : length(x), ]
    singleton_responses <- as.list(predicted_nn[, 1])
    singleton_parents <- lapply(seq_len(nrow(predicted_nn)), function(i) {
      output <- predicted_nn[i, -1]
      output <- output[!is.na(output)]
      output <- output[order(x[output])]
      output
    })
  }

  block_responses <- c(sequential_block_responses, singleton_responses)
  block_parents <- c(sequential_block_parents, singleton_parents)
  block_include_noise <- lapply(seq_along(block_parents), function(i) {
    c(block_parents[[i]], block_responses[[i]]) <= n_observed
  })

  observed_indices <- 1 : (d * n_observed)
  predicted_indices <- (d * n_observed + 1) : (d * length(x))
  if (has_cache) {
    global_hash <- rlang::hash(list(
      fit = fit,
      block_id = .get_block_ids(
        x,
        block_parents,
        block_responses,
        block_include_noise,
        cache_resolution
      ),
      observed_indices = observed_indices,
      predicted_indices = predicted_indices
    ))
    global_hash_U <- sprintf('u000%s', global_hash)
  }

  if (!has_cache || !cache$exists(global_hash_U)) {
    U_parts <- .make_U_parts(
      fit,
      x,
      d,
      block_parents,
      block_responses,
      block_include_noise,
      local_cache,
      cache_resolution
    )
    U <- Matrix::sparseMatrix(
      i = do.call(c, lapply(U_parts, getElement, 'row_indices')),
      j = do.call(c, lapply(U_parts, getElement, 'col_indices')),
      x = do.call(c, lapply(U_parts, getElement, 'entries')),
      triangular = TRUE
    )
    U_oo <- U[observed_indices, observed_indices]
    U_pp <- U[predicted_indices, predicted_indices]
    U_pa <- U[predicted_indices, ]
    U_oa <- U[observed_indices, ]

    if (has_cache) {
      cache$set(global_hash_U, list(
        U_oo = U_oo,
        U_pp = U_pp,
        U_pa = U_pa,
        U_oa = U_oa
      ))
    }
  } else if (has_cache) {
    logger::log_trace('Loading U matrix from cache')
    cache_entry_U <- cache$get(global_hash_U)
    U_oo <- cache_entry_U$U_oo
    U_pp <- cache_entry_U$U_pp
    U_pa <- cache_entry_U$U_pa
    U_oa <- cache_entry_U$U_oa
  }

  logger::log_trace('Calculating beta_hat')
  X_o <- observed_model_data_long$X
  X_p <- predicted_model_data_long$X
  y_o <- observed_model_data_long$y

  Xt_Q_X <- crossprod(as.matrix(Matrix::crossprod(U_oo, X_o)))
  Q_y <- as.vector(U_oo %*% Matrix::crossprod(U_oo, y_o))
  beta_hat <- as.vector(solve(Xt_Q_X, crossprod(X_o, Q_y)))

  logger::log_trace('Calculating fitted and mean')
  y_o_tilde <- as.vector(y_o - X_o %*% beta_hat)
  y_predicted_fitted <- as.vector(X_p %*% beta_hat)
  y_predicted_mean <- y_predicted_fitted - as.vector(
    Matrix::solve(
      Matrix::t(U_pp),
      Matrix::solve(
        U_pp,
        as.vector(U_pa %*% Matrix::crossprod(U_oa, y_o_tilde))
      )
    )
  )

  if (include_posterior_covariance) {
    if (has_cache) {
      global_hash_Sigma_y_prediction <- sprintf('s000%s', global_hash)
    }

    if (!has_cache || !cache$exists(global_hash_Sigma_y_prediction)) {
      logger::log_trace('Calculating prediction covariance')
      Sigma_y_predicted <- as.matrix(solve(as.matrix(Matrix::tcrossprod(U_pp))))

      if (has_cache) {
        cache$set(global_hash_Sigma_y_prediction, Sigma_y_predicted)
      }
    } else if (has_cache) {
      logger::log_trace('Loading prediction covariance from cache')
      Sigma_y_predicted <- cache$get(global_hash_Sigma_y_prediction)
    }
  }

  if (n_samples > 0) {
    logger::log_trace('Sampling from prediction distribution')
    y_tilde <- as.matrix(Matrix::solve(
      Matrix::t(U_pp),
      matrix(rnorm(n_samples * n_predicted * d), ncol = n_samples)
    ))
    y_samples <- y_predicted_mean + y_tilde
  }

  list(
    fitted = y_predicted_fitted,
    mean = y_predicted_mean,
    covariance = if (include_posterior_covariance) {
      Sigma_y_predicted
    } else {
      NULL
    },
    samples = if (n_samples > 0) y_samples else NULL
  )
}

.ptide_predict_general <- function(
  fit,
  observed_model_data,
  predicted_model_data,
  observed_model_data_long,
  predicted_model_data_long,
  include_posterior_covariance,
  n_samples,
  vecchia_configuration,
  cache,
  cache_resolution
) {
  n_parents <- vecchia_configuration$n_parents
  n_responses <- vecchia_configuration$n_responses
  d <- observed_model_data$d
  n_observed <- length(observed_model_data$x)
  n_predicted <- length(predicted_model_data$x)
  has_cache <- !missing(cache)

  local_cache <- if (has_cache) {
    cache
  } else {
    cachem::cache_mem(evict = 'fifo')
  }

  if (has_cache) {
    min_x <- min(c(min(observed_model_data$x), min(predicted_model_data$x)))
    global_hash <- rlang::hash(list(
      fit = fit,
      vecchia_configuration = vecchia_configuration,
      observed_x = round(cache_resolution * (observed_model_data$x - min_x)),
      predicted_x = round(cache_resolution * (predicted_model_data$x - min_x))
    ))
    global_hash_U <- sprintf('ug000%s', global_hash)
  }

  x_unique <- sort(unique(c(observed_model_data$x, predicted_model_data$x)))

  if (!has_cache || !cache$exists(global_hash_U)) {
    n_latent_blocks <- ceiling((length(x_unique) - n_parents) / n_responses)
    latent_block_responses <- c(
      list(1L : (n_parents + n_responses)),
      lapply(2L : n_latent_blocks, function(i) {
        ((i - 1L) * n_responses + n_parents + 1L) : min(
          i * n_responses + n_parents,
          length(x_unique)
        )
      })
    )
    latent_block_parents <- c(
      list(integer(0)),
      lapply(2L : n_latent_blocks, function(i) {
        ((i - 1L) * n_responses + 1L) : ((i - 1L) * n_responses + n_parents)
      })
    )
    latent_block_include_noise <- lapply(seq_len(n_latent_blocks), function(i) {
      rep(
        FALSE,
        length(latent_block_parents[[i]]) + length(latent_block_responses[[i]])
      )
    })

    latent_U_parts <- .make_U_parts(
      fit,
      x_unique,
      d,
      latent_block_parents,
      latent_block_responses,
      latent_block_include_noise,
      local_cache,
      cache_resolution
    )

    observed_Sigma_base_array <- ptide_covariance(
      c(0, 0),
      fit,
      c(0, 1)
    )
    observed_Sigma_base_matrix <- matrix(
      aperm(observed_Sigma_base_array, c(3, 1, 4, 2)),
      nrow = prod(dim(observed_Sigma_base_array)[c(1, 3)])
    )
    observed_R <- chol(observed_Sigma_base_matrix)
    observed_B <- rbind(
      matrix(0, nrow = d, ncol = d),
      diag(d)
    )
    observed_values <- as.vector(backsolve(observed_R, observed_B))
    observed_values_non_zero <- observed_values != 0
    observed_U_parts <- lapply(seq_len(n_observed), function(i) {
      parent_index <- match(observed_model_data$x[i], x_unique)
      response_index <- length(x_unique) + i
      indices_j <- c(parent_index, response_index)
      row_indices_j <- (
        rep(d * indices_j, each = d)
        + rep((-d + 1L) : 0L, length(indices_j))
      )
      col_indices_j <- (
        rep(d * response_index, each = d)
        + rep((-d + 1L) : 0L, 1)
      )
      list(
        row_indices = rep(
          row_indices_j,
          length(col_indices_j)
        )[observed_values_non_zero],
        col_indices = rep(
          col_indices_j,
          each = length(row_indices_j)
        )[observed_values_non_zero],
        entries = observed_values[observed_values_non_zero]
      )
    })

    join_parts <- function(name){
      c(
        do.call(c, lapply(latent_U_parts, getElement, name)),
        do.call(c, lapply(observed_U_parts, getElement, name))
      )
    }

    U <- Matrix::sparseMatrix(
      i = join_parts('row_indices'),
      j = join_parts('col_indices'),
      x = join_parts('entries'),
      triangular = TRUE
    )
    U_l <- U[seq_len(d * length(x_unique)), ]
    U_o <- U[(d * length(x_unique) + 1) : nrow(U), ]

    logger::log_trace('Forming W')
    W <- Matrix::tcrossprod(U_l)

    logger::log_trace('Calculating V')
    rev_i <- rev(seq_len(nrow(W)))
    V_a <- Matrix::t(Matrix::chol(W[rev_i, rev_i]))
    V <- as(V_a[rev_i, rev_i], 'triangularMatrix')

    if (has_cache) {
      cache$set(global_hash_U, list(
        U_l = U_l,
        U_o = U_o,
        V = V
      ))
    }
  } else if (has_cache) {
    logger::log_trace('Loading U matrix from cache')
    cache_entry <- cache$get(global_hash_U)
    U_l <- cache_entry$U_l
    U_o <- cache_entry$U_o
    V <- cache_entry$V
  }

  X_o <- observed_model_data_long$X
  z_o <- observed_model_data_long$y
  X_p <- predicted_model_data_long$X

  X_latent <- matrix(NA, nrow = d * length(x_unique), ncol = ncol(X_o))
  latent_in_observed <- match(x_unique, observed_model_data$x)
  latent_in_predicted <- match(x_unique, predicted_model_data$x)
  predicted_in_latent <- match(predicted_model_data$x, x_unique)

  latent_in_observed_long <- (
    rep(d * latent_in_observed, each = d)
    + rep((-d + 1L) : 0L, length(latent_in_observed))
  )
  latent_in_predicted_long <- (
    rep(d * latent_in_predicted, each = d)
    + rep((-d + 1L) : 0L, length(latent_in_predicted))
  )
  predicted_in_latent_long <- (
    rep(d * predicted_in_latent, each = d)
    + rep((-d + 1L) : 0L, length(predicted_in_latent))
  )

  X_latent[!is.na(latent_in_observed_long), ] <- X_o[
    latent_in_observed_long[!is.na(latent_in_observed_long)],
  ]
  X_latent[!is.na(latent_in_predicted_long), ] <- X_p[latent_in_predicted_long[
    !is.na(latent_in_predicted_long)],
  ]
  stopifnot(!anyNA(X_latent))

  if (has_cache) {
    global_hash_Q_X <- sprintf('qg000%s', global_hash)
  }
  if (!has_cache || !cache$exists(global_hash_Q_X)) {
    logger::log_trace('Calculating Q_X')
    Q_X <- (
      as.matrix(U_o %*% Matrix::crossprod(U_o, X_o))
      - as.matrix(
        U_o %*% Matrix::crossprod(U_l, Matrix::solve(
          Matrix::t(V), Matrix::solve(V, U_l %*% Matrix::crossprod(U_o, X_o))
        ))
      )
    )
    if (has_cache) {
      cache$set(global_hash_Q_X, Q_X)
    }
  } else if (has_cache) {
    Q_X <- cache$get(global_hash_Q_X)
  }

  logger::log_trace('Calculating Q_z')
  Q_z <- (
    as.vector(U_o %*% Matrix::crossprod(U_o, z_o))
    - as.vector(
      U_o %*% Matrix::crossprod(U_l, Matrix::solve(
        Matrix::t(V), Matrix::solve(V, U_l %*% Matrix::crossprod(U_o, z_o))
      ))
    )
  )

  logger::log_trace('Calculating beta_hat')
  Xt_Q_X <- crossprod(X_o, Q_X)
  beta_hat <- as.vector(solve(Xt_Q_X, crossprod(X_o, Q_z)))

  logger::log_trace('Calculating y_latent_mean')
  z_o_tilde <- as.vector(z_o - X_o %*% beta_hat)
  y_latent_fitted <- as.vector(X_latent %*% beta_hat)
  y_latent_mean <- y_latent_fitted - as.vector(
    Matrix::solve(Matrix::t(V), Matrix::solve(
      V,
      U_l %*% Matrix::crossprod(U_o, z_o_tilde)
    ))
  )

  y_predicted_fitted <- y_latent_fitted[predicted_in_latent_long]
  y_predicted_mean <- y_latent_mean[predicted_in_latent_long]

  if (include_posterior_covariance || n_samples > 0) {
    if (has_cache) {
      global_hash_Sigma_y_prediction <- sprintf('sg000%s', global_hash)
    }

    if (!has_cache || !cache$exists(global_hash_Sigma_y_prediction)) {
      logger::log_trace('Calculating prediction covariance')
      H_t_dense <- matrix(
        0,
        nrow = d * length(x_unique),
        ncol = d * n_predicted
      )
      H_t_dense[cbind(
        predicted_in_latent_long,
        seq_len(d * n_predicted)
      )] <- 1
      Sigma_y_predicted <- crossprod(as.matrix(Matrix::solve(V, H_t_dense)))

      if (has_cache) {
        cache$set(global_hash_Sigma_y_prediction, Sigma_y_predicted)
      }
    } else if (has_cache) {
      Sigma_y_predicted <- cache$get(global_hash_Sigma_y_prediction)
    }
  }

  if (n_samples > 0) {
    logger::log_trace('Drawing samples')
    y_tilde <- as.matrix(crossprod(
      chol(Sigma_y_predicted),
      matrix(rnorm(n_samples * d * n_predicted), ncol = n_samples)
    ))
    y_samples <- y_predicted_mean + y_tilde
  }

  list(
    fitted = y_predicted_fitted,
    mean = y_predicted_mean,
    covariance = if (include_posterior_covariance) {
      Sigma_y_predicted
    } else {
      NULL
    },
    samples = if (n_samples > 0) y_samples else NULL
  )
}

.get_block_ids <- function(
  x,
  block_parents,
  block_responses,
  block_include_noise,
  cache_resolution
) {
  lapply(seq_along(block_parents), function(i) {
    block_x <- x[c(block_parents[[i]], block_responses[[i]])]
    list(
      x_rounded = round(cache_resolution * (block_x - min(block_x))),
      include_noise = block_include_noise[[i]],
      n_parents = length(block_parents[[i]]),
      n_responses = length(block_responses[[i]])
    )
  })
}

.make_U_parts <- function(
  fit,
  x,
  d,
  block_parents,
  block_responses,
  block_include_noise,
  cache,
  cache_resolution
) {
  block_id <- .get_block_ids(
    x,
    block_parents,
    block_responses,
    block_include_noise,
    cache_resolution
  )
  block_id_hash <- sapply(block_id, rlang::hash)
  block_id_hash_table_full <- table(block_id_hash)
  block_kernel_group <- match(block_id_hash, names(block_id_hash_table_full))
  n_kernel_groups <- max(block_kernel_group)

  logger::log_trace(
    'Calculating U parts using {length(block_id)} blocks across',
    ' {n_kernel_groups} groups'
  )

  lapply(seq_len(n_kernel_groups), function(i) {
    blocks <- which(block_kernel_group == i)
    x_i <- x[c(block_parents[[blocks[1]]], block_responses[[blocks[1]]])]
    n_parents_i <- length(block_parents[[blocks[1]]])
    n_responses_i <- length(block_responses[[blocks[1]]])
    include_noise_i <- block_include_noise[[blocks[1]]]
    id_R <- list(
      fit = fit,
      x_rounded = round(cache_resolution * (x_i - min(x_i))),
      include_noise = include_noise_i
    )
    id_R_hash <- sprintf('r000%s', rlang::hash(id_R))
    if (!cache$exists(id_R_hash)) {
      Sigma_array <- ptide_covariance(
        x_i,
        fit,
        include_noise_i
      )
      Sigma_matrix <- matrix(
        aperm(Sigma_array, c(3, 1, 4, 2)),
        nrow = prod(dim(Sigma_array)[c(1, 3)])
      )
      R <- chol(Sigma_matrix)
      cache$set(id_R_hash, R)
    } else {
      R <- cache$get(id_R_hash)
    }
    B <- rbind(
      matrix(0, nrow = d * n_parents_i, ncol = d * n_responses_i),
      diag(d * n_responses_i)
    )
    values <- as.vector(backsolve(R, B))
    values_non_zero <- values != 0

    parts <- lapply(blocks, function(j) {
      indices_j <- c(block_parents[[j]], block_responses[[j]])

      row_indices_j <- (
        rep(d * indices_j, each = d)
        + rep((-d + 1L) : 0L, length(indices_j))
      )
      col_indices_j <- (
        rep(d * block_responses[[j]], each = d)
        + rep((-d + 1L) : 0L, length(block_responses[[j]]))
      )
      list(
        row_indices = rep(
          row_indices_j,
          length(col_indices_j)
        )[values_non_zero],
        col_indices = rep(
          col_indices_j,
          each = length(row_indices_j)
        )[values_non_zero],
        entries = values[values_non_zero]
      )
    })
    list(
      row_indices = do.call(c, lapply(parts, getElement, 'row_indices')),
      col_indices = do.call(c, lapply(parts, getElement, 'col_indices')),
      entries = do.call(c, lapply(parts, getElement, 'entries'))
    )
  })
}
