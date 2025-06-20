#' Perform time series cross-validation for PTide
#'
#' Perform time series cross-validatiom for a PTide fit. Given a dataset,
#' this method will perform a series of predictions by advancing the
#' observed data by \code{skip} observation at a time.
#'
#' @param fit A PTide fit object as returned by \code{\link{ptide_fit}}.
#' @param data A data frame containing the data to be used for the validation.
#' @param start_index The row number of the first observation from which to
#' start the validation.
#' @param horizons A vector of time horizons for which to make predictions.
#' Should be a period of time created using \code{\link{lubridate::hours}} or
#' similar.
#' @param skip The number of observations to skip between each prediction.
#' @param last_n The number of observations to condition on for the prediction.
#' If \code{Inf}, conditions on all previous observations.
#' @param include_sd Whether to include the standard deviation of the forecast
#' in the output.
#' @param include_marginal_correlation Whether to include the marginal
#' correlation of the forecast in the output.
#' @param cores The number of cores to use for parallelisation. If \code{1},
#' parallelisation is disabled.
#' @param cache A cache object satisfying the semantics of the \code{cachem}
#' package, such as \code{\link[cachem]{cache_mem}}. This is passed to
#' \code{\link{ptide_predict}}.
#' @param ... Additional arguments to be passed to \code{\link{ptide_predict}}.
#'
#' @return A data frame containing the predictions alongside the input data.
#' @export
run_validation <- function(
  fit,
  data,
  start_index,
  horizons,
  skip = 1L,
  last_n = Inf,
  include_sd = FALSE,
  include_marginal_correlation = FALSE,
  cores = 1L,
  cache = cachem::cache_mem(),
  ...
) {
  stopifnot(!is.unsorted(data[[fit$model$time_field]]))

  indices <- seq(start_index, nrow(data), by = skip)

  do_jobs <- if (cores == 1L) {
    function(...) lapply(...)
  } else {
    function(...) pbmcapply::pbmclapply(
      ...,
      mc.cores = cores,
      ignore.interactive = TRUE
    )
  }

  print_error <- function(fn) {
    function(...) {
      tryCatch(
        fn(...),
        error = function(e) {
          print(e)
          stop(e)
        }
      )
    }
  }

  input_hours <- .model_hours(fit$model, data)

  parts <- do_jobs(indices, print_error(function(index) {
    training_start_index <- max(1L, index - last_n + 1L)

    time_issued <- data[[fit$model$time_field]][index]

    forecast_times <- time_issued + horizons
    forecast_hours <- .model_hours(fit$model, time = forecast_times)

    forecast_data <- data[input_hours %in% forecast_hours, ]
    if (nrow(forecast_data) == 0) return(NULL)

    forecast <- ptide_predict(
      fit,
      forecast_data,
      data[training_start_index : index, ],
      include_sd = include_sd,
      include_marginal_covariance = include_marginal_correlation,
      cache = cache,
      ...
    )
    forecast_mean <- forecast$mean
    colnames(forecast_mean) <- paste0(colnames(forecast_mean), '_mean')

    forecast_data$time_issued <- time_issued
    colnames(forecast_data)[colnames(forecast_data) == fit$model$time_field] <- 'time_forecast'

    output <- cbind(forecast_data, forecast_mean)

    if (include_sd) {
      forecast_sd <- forecast$sd
      colnames(forecast_sd) <- paste0(colnames(forecast_sd), '_sd')
      output <- cbind(output, forecast_sd)
    }

    if (include_marginal_correlation) {
      n <- dim(forecast$marginal_covariance)[1]
      d <- dim(forecast$marginal_covariance)[2]
      if (d > 1) {
        forecast_correlation <- matrix(0, nrow = n, ncol = d * (d - 1) / 2)
        k <- 1
        for (i in 1 : (d - 1)) {
          for (j in (i + 1) : d) {
            forecast_correlation[, k] <- forecast$marginal_covariance[, i, j] / sqrt(
              forecast$marginal_covariance[, i, i] * forecast$marginal_covariance[, j, j]
            )
            colnames(forecast_correlation)[k] <- paste0('cor_', paste0(
              dimnames(forecast$marginal_covariance)[[2]][c(i, j)],
              collapse = '_'
            ))
            k <- k + 1
          }
        }

        output <- cbind(output, forecast_correlation)
      }
    }

    output
  }))

  any_error <- FALSE
  for (i in seq_along(parts)) {
    if (inherits(parts[[i]], 'try-error')) {
      cat('part', i, 'had error:\n')
      print(parts[[i]])
      any_error <- TRUE
    }
  }
  if (any_error) {
    stop('Error computing validation')
  }

  dplyr::bind_rows(parts)
}
