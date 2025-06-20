#' Specify a PTide model
#'
#' Describe a PTide model, including response variable, external covariates,
#' tidal constituents, and covariance functions. This function allows the user
#' to specify the structure of a PTide model, including its formula, start time,
#' tidal constituents, mean function, covariance functions, and the time field.
#'
#' @param formula An object of class \code{\link[stats]{formula}} specifying the
#' model to be fitted. This formula includes the response variable or variables
#' and potentially external covariates. If more than one response variable is
#' specified (using \link{cbind}, the model is a multivariate model.
#' @param start_time A POSIXct object representing the start time of the time
#' series. This is used to align the model with the time field of the data. If
#' not provided, a default of '2020-01-01 00:00:00' is used.
#' @param tidal_constituents A character vector specifying the tidal
#' constituents to be included in the model. Default constituents are 'M2' and
#' 'S2'.
#' @param mean_function A list specifying the prior for the coefficients of the
#' mean function of the model. The list should contain two elements: 'mean'
#' (the prior mean), either a scalar or a vector, and 'precision' (the prior
# precision, the reciprocal of the variance), either a scalar or a matrix.
#' Default values are 0 for mean and 1/100 for precision.
#' @param covariance_functions A list of functions specifying the covariance
#' structure of the model. Each function in the list should be one of the
#' covariance functions defined in the package (like
#' \code{\link{ptide_exponential_covariance_function}}). By default,
#' exponential, quasi-periodic, and white noise covariance functions are
#' included.
#' @param time_field A string specifying the name of the field in the data that
#" represents time.
#'
#' @return An object of class 'ptide_model' that represents the specified PTide
#' model.
#'
#' @examples
#' # Bivariate response
#' ptide_model(
#'   formula = cbind(easting, northing) ~ 1,
#'   start_time = lubridate::ymd_hms('2023-01-01 00:00:00'),
#'   tidal_constituents = c('M2', 'S2', 'K1')
#' )
#'
#' @export
ptide_model <- function(
  formula,
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
  time_field = 'time'
) {
  if (missing(start_time)) {
    start_time <- lubridate::ymd_hms('2020-01-01 00:00:00')
  }
  stopifnot(length(covariance_functions) > 0)
  structure(list(
    formula = formula,
    start_time = start_time,
    tidal_constituents = tidal_constituents,
    mean_function = mean_function,
    covariance_functions = covariance_functions,
    time_field = time_field
  ), class = 'ptide_model')
}

.model_hours <- function(model, data, time = data[[model$time_field]]) {
  as.double(difftime(
    time,
    model$start_time,
    units = 'hours'
  ))
}
