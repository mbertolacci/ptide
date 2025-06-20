#' Specify a PTide model and fit it to data
#'
#' Specify and fit a PTide model to data.
#'
#' @param formula A formula specifying the model; see \code{\link{ptide_model}}
#' for more information.
#' @param data A data frame containing the data to fit the model to.
#' @param start_time A \code{POSIXct} object specifying the start time of the
#' model; see \code{\link{ptide_model}} for more information.
#' @param tidal_constituents A vector of strings specifying the tidal
#' constituents to include in the model; see \code{\link{ptide_model}} for more
#' information.
#' @param mean_function A list specifying the prior for the mean function.
#' See \code{\link{ptide_model}} for more information.
#' @param covariance_functions A list of functions specifying the covariance
#' structure of the model. See \code{\link{ptide_model}} for more information.
#' @param time_field A string specifying the name of the field in the data that
#' represents time.
#' @param ... Additional arguments passed to \code{\link{ptide_fit}}.
#' @export
ptide <- function(
  formula,
  data,
  start_time,
  tidal_constituents = c('M2', 'S2'),
  mean_function = list(
    mean = 0,
    precision = 1 / 100
  ),
  covariance_functions = list(
    ptide_exponential_covariance_function(),
    ptide_quasi_periodic_covariance_function(),
    ptide_white_noise_covariance_function()
  ),
  time_field = 'time',
  ...
) {
  model <- ptide_model(
    formula,
    start_time,
    tidal_constituents,
    mean_function,
    covariance_functions,
    time_field
  )

  ptide_fit(
    model,
    data,
    ...
  )
}
