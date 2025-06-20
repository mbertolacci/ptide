#' Fit a PTide model to data
#'
#' Fits a PTide model specified with \code{\link{ptide_model}} to data.
#'
#' The computational time required for the fit depends on the number of
#' observations, the choice of Vecchia approximation (see the
#' \code{vecchia_configuration} parameter below), and whether the observed
#' times are regular. The last point is particularly important: if the observed
#" times are regular, many matrix operations can be precomputed and the
#' optimisation routine will run much faster.
#'
#' @param model A PTide model specified with \code{\link{ptide_model}}.
#' @param data A data frame containing the data to fit the model to.
#' @param vecchia_configuration A list containing the configuration for the
#' Vecchia approximation. The list must contain two elements: \code{n_parents}
#' (number of parents) and \code{n_responses} (number of responses). For the
#' approximation, the time series is split into blocks of length
#' \code{n_responses}, and each block conditions on the previous
#' \code{n_parents} time points.
#' @param max_attempts Maximum number of optimisation attempts.
#' @param best_of Number of times to run the optimisation algorithm. The
#' algorithm is run repeatedly with different starting values in order to help
#' find the global posterior mode.
#' @param grain_size Used to control the parallelisation. If set to 'auto',
#' this will attempted to maximise the use of threads
#' @param threads Number of threads to use for parallelisation, taken from
#' `getOption('ptide.threads')` by default, which is itself set to -1 by
#' default. The special value -1 picks a number of threads based on the number
#' of cores in the system).
#' @param x Object to be printed in \code{print.ptide_fit}.
#' @param ... Additional parameters to be passed to
#' \code{\link{stan_optimizing_best_of}}
#' @return A PTide fit object (a list) with the following entries:
#' \itemize{
#'   \item \code{model}: The PTide model object.
#'   \item \code{covariance_functions}: A list of covariance function objects,
#'     one for each covariance function in the model, with the fitted
#'     parameters.
#'   \item \code{beta_hat}: The fitted mean function parameters.
#'   \item \code{optimization_result}: The raw output of the optimisation
#'     procedure.
#' }
#' @export
ptide_fit <- function(
  model,
  data,
  vecchia_configuration = list(n_parents = 256L, n_responses = 256L),
  max_attempts = 20,
  best_of = 10,
  grain_size = 'auto',
  threads = getOption('ptide.threads', -1L),
  ...
) {
  stopifnot(!is.unsorted(data[[model$time_field]]))

  .local_tbb_threads(threads)

  stan_data <- ptide_stan_data(model, data, vecchia_configuration)
  stan_data$use_parallel <- threads == -1L || threads > 1L
  stan_data$strategy <- 0L

  if (grain_size == 'auto') {
    if (stan_data$use_parallel) {
      n_threads <- .get_tbb_threads()
      grain_size <- ceiling(stan_data$n_blocks / n_threads)
      logger::log_trace(
        'Number of threads = {n_threads}, ',
        'number of blocks = {stan_data$n_blocks}, ',
        'chosen grain size = {grain_size}'
      )
    } else {
      grain_size <- 1L
    }
  }
  stan_data$grain_size <- grain_size

  optimization_result <- stan_optimizing_best_of(
    object = ptide_stan_model(),
    data = stan_data,
    .max_attempts = max_attempts,
    .best_of = best_of,
    .method = 'stan_nlminb',
    ...
  )

  # Unpack optimised kernel hyperparameters into output object
  covariance_functions <- do.call(c, lapply(
    .covariance_function_types,
    function(type) {
      output <- Filter(function(x) x$type == type, model$covariance_functions)
      if (length(output) == 0) return(NULL)
      parameter_names <- c(
        'sigma_squared',
        'L',
        .covariance_function_parameter_names[[type]]
      )
      for (i in seq_along(output)) {
        output[[i]] <- output[[i]][c(
          'type',
          .covariance_function_configuration_names[[type]]
        )]
        for (parameter_name in parameter_names) {
          output[[i]][[parameter_name]] <- .index_array(
            optimization_result$par[[sprintf('%s_%s', type, parameter_name)]],
            1,
            i,
            drop = TRUE
          )
        }
      }
      output
    }
  ))

  beta_hat_matrix <- matrix(
    optimization_result$par$beta_hat,
    nrow = ncol(stan_data$X)
  )
  colnames(beta_hat_matrix) <- colnames(stan_data$y)

  structure(list(
    model = model,
    covariance_functions = covariance_functions,
    beta_hat = beta_hat_matrix,
    optimization_result = optimization_result
  ), class = 'ptide_fit')
}

#' @describeIn ptide_fit Print a fit to screen
#' @export
print.ptide_fit <- function(x, ...) {
  printf <- function(...) cat(sprintf(...))

  cat('-- Formula\n')
  print(x$model$formula)

  cat('-- Tidal constituents and amplitudes (period in hours)\n')
  n_constituents <- length(x$model$tidal_constituents)
  d <- ncol(x$beta_hat)
  harmonic_beta_hat <- tail(x$beta_hat, 2L * n_constituents)
  amplitudes <- matrix(0, nrow = n_constituents, ncol = d)
  for (i in seq_len(n_constituents)) {
    amplitudes[i, ] <- sqrt(
      harmonic_beta_hat[2L * i - 1L, ] ^ 2
      + harmonic_beta_hat[2L * i, ] ^ 2
    )
  }
  colnames(amplitudes) <- sprintf(
    '%s_amplitude',
    colnames(x$beta_hat)
  )
  print(
    ptide::tidal_constituents %>%
      dplyr::filter(name %in% x$model$tidal_constituents) %>%
      dplyr::select(name, freq) %>%
      dplyr::mutate(period = 1 / freq) %>%
      cbind(amplitudes)
  )

  cat('-- Covariance functions\n')
  for (covariance_function in x$covariance_functions) {
    printf('- Type: %s\n', covariance_function$type)
    variable_names <- .covariance_function_configuration_names[[
      covariance_function$type
    ]]
    parameter_names <- .covariance_function_parameter_names[[
      covariance_function$type
    ]]
    if (length(variable_names) > 0) {
      cat('  Configuration:\n')
      for (variable_name in variable_names) {
        printf(
          '    %s = %f\n',
          variable_name, covariance_function[[variable_name]]
        )
      }
    }
    cat('  Parameters:\n')
    for (parameter_name in parameter_names) {
      printf(
        '    %s = %f\n',
        parameter_name,
        covariance_function[[parameter_name]]
      )
    }
    cat('    Sigma =\n')
    print(tcrossprod(
      diag(sqrt(covariance_function$sigma_squared)) %*% covariance_function$L
    ))
  }
}

.ptide_model_data <- function(
  model,
  data,
  include_response = TRUE
) {
  terms <- terms(model$formula)
  if (!include_response) {
    terms <- delete.response(terms)
  }

  hours <- .model_hours(model, data)
  model_frame <- model.frame(terms, data, na.action = 'na.fail')

  output <- list(
    X = cbind(
      model.matrix(terms, model_frame),
      as.matrix(tidal_harmonic_design_matrix(
        hours = hours,
        constituent_names = model$tidal_constituents
      ))
    ),
    x = hours
  )
  if (include_response) {
    output$y <- as.matrix(model.response(model_frame))
    output$d <- ncol(output$y)
  } else {
    output$d <- length(all.vars(model$formula)) - length(all.vars(terms))
  }
  output
}

#' @export
ptide_stan_data <- function(model, data, vecchia_configuration, ...) {
  output <- .ptide_model_data(model, data, ...)

  has_empirical_variance <- any(sapply(
    model$covariance_functions,
    getElement,
    'sigma_squared_a'
  ) == 'empirical')
  if (has_empirical_variance) {
    variances <- sapply(seq_len(output$d), function(i) {
      X <- output$X
      y <- output$y[, i]
      var(y - X %*% solve(crossprod(X), crossprod(X, y)))
    })
    variance_parameters <- inverse_gamma_quantile_prior(
      min(variances) / 10,
      max(variances) * 10,
      0.001,
      0.999
    )
  }

  covariance_function_data <- do.call(c, lapply(
    .covariance_function_types,
    function(type) {
      variable_names <- c(
        .covariance_function_configuration_names[[type]],
        .covariance_function_hyperparameter_names[[type]]
      )
      covariance_functions_of_type <- Filter(
        function(x) x$type == type,
        model$covariance_functions
      )

      output <- list()
      output[[sprintf('n_%s_kernels', type)]] <- length(
        covariance_functions_of_type
      )
      for (variable_name in variable_names) {
        output_name <- sprintf('%s_%s', type, variable_name)
        if (length(covariance_functions_of_type) == 0) {
          output[[output_name]] <- array(0, dim = c(0))
          next
        }
        output[[output_name]] <- array(if (
          has_empirical_variance
          && startsWith(variable_name, 'sigma_squared_')
        ) {
          sapply(covariance_functions_of_type, function(covariance_function) {
            if (covariance_function$sigma_squared_a == 'empirical') {
              if (variable_name == 'sigma_squared_a') {
                variance_parameters$shape
              } else {
                variance_parameters$rate
              }
            } else {
              covariance_function[[variable_name]]
            }
          })
        } else {
          sapply(covariance_functions_of_type, getElement, variable_name)
        })
      }
      output
    }
  ))

  beta_mean <- if (length(model$mean_function$mean) > 1) {
    model$mean_function$mean
  } else {
    rep(model$mean_function$mean, ncol(output$X) * output$d)
  }
  beta_precision <- if (is.matrix(model$mean_function$precision)) {
    model$mean_function$precision
  } else {
    diag(rep(model$mean_function$precision, ncol(output$X) * output$d))
  }

  c(
    output,
    list(
      n = nrow(output$y),
      p = ncol(output$X),
      beta_prior_mean = beta_mean,
      beta_prior_precision = beta_precision
    ),
    covariance_function_data,
    .vecchia_blocks(output$x, vecchia_configuration)
  )
}

.vecchia_blocks <- function(x, vecchia_configuration) {
  n_parents <- vecchia_configuration$n_parents
  n_responses <- vecchia_configuration$n_responses
  n <- length(x)
  block_size <- n_parents + n_responses
  # First block has no parents
  n_blocks <- ceiling((n - block_size) / n_responses) + 1L
  block_indices <- seq_len(block_size)
  block_last_index <- block_size
  block_n_responses <- block_size
  for (i in seq_len(n_blocks - 1L)) {
    start_response_i <- block_size + (i - 1L) * n_responses + 1L
    n_responses_i <- min(n_responses, n - start_response_i + 1L)
    block_indices <- c(
      block_indices,
      (start_response_i - n_parents) : (start_response_i + n_responses_i - 1L)
    )
    block_last_index <- c(block_last_index, length(block_indices))
    block_n_responses <- c(block_n_responses, n_responses_i)
  }

  block_x <- lapply(seq_len(n_blocks), function(i) {
    if (i == 1) {
      start_index <- 1L
    } else {
      start_index <- block_last_index[i - 1L] + 1L
    }
    end_index <- block_last_index[i]
    output <- x[block_indices[start_index : end_index]]
    round(60 * output - 60 * min(output))
  })

  block_x_hash <- sapply(block_x, rlang::hash)
  block_x_hash_table_full <- table(block_x_hash)
  block_kernel_group <- match(block_x_hash, names(block_x_hash_table_full))

  list(
    n_indices = length(block_indices),
    n_blocks = n_blocks,
    block_indices = array(as.integer(block_indices)),
    block_last_index = array(as.integer(block_last_index)),
    block_n_responses = array(as.integer(block_n_responses)),
    block_kernel_group = array(as.integer(block_kernel_group))
  )
}

.lengthen_tidal_model_data <- function(data, stack = TRUE) {
  output <- list(
    x = data$x,
    d = data$d
  )
  if (stack) {
    output$X <- kronecker(diag(data$d), data$X)
    if (!is.null(data$y)) {
      output$y <- as.vector(data$y)
    }
  } else {
    output$X <- kronecker(data$X, diag(data$d))
    if (!is.null(data$y)) {
      output$y <- as.vector(t(data$y))
    }
  }
  output
}

#' @export
ptide_stan_model <- function() {
  stanmodels$ptide
}
