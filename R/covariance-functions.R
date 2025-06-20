#' Covariance functions
#'
#' These functions are used to specify covariance functions for use with
#' \code{\link{ptide}} and \code{\link{ptide_model}}.
#'
#' The covariance functions in this package have a separable structure, with
#' \eqn{K(t, t') = k(t, t')S} a \eqn{d \times d} matrix, where d is the
#' dimension of the process, \eqn{k(t, t')} is a univariate covariance function,
#' and S is a \eqn{d \times d} covariance matrix. Internally, S is parameterised
#' as diag(sigma) %*% C %*% diag(sigma), where sigma is a \eqn{d \times 1}
#' vector of standard deviations, and C is a correlation matrix. The univariate
#' covariate functions available are described below. All length scale
#' parameters are given inverse gamma priors, sigma^2 is given an inverse gamma
#' prior, and C is given the LKJ prior.
#' @section Exponential:
#' The exponential covariance function
#' (\code{ptide_exponential_covariance_function}) has the expression
#' \deqn{k(t, t') = \exp\left( -\frac{|t - t'|}{\ell / d} \right),}
#' where \eqn{\ell} is the length scale parameter, and d is a
#' user-specified length dilation (called \code{length_dilation} in the
#' parameters).
#' @section Matern 3/2 and 5/2:
#' The Matern 3/2 and 5/2 covariance functions
#' (\code{ptide_matern32_covariance_function} and
#' \code{ptide_matern52_covariance_function}) has the expression
#' \deqn{
#'   k(t, t') = \frac{
#'     2^{1 - \nu}}{\Gamma(\nu)
#'   }(\sqrt{2\nu}|t - t'| / \ell)^\nu K_\nu(\sqrt{2 \nu} |t - t'| / \ell),}
#' where \eqn{\nu = 3/2} or \eqn{5 / 2}, \eqn{\ell} is the length scale
#' parameter, and d is a user-specified length dilation.
#' @section Squared exponential:
#' The squared exponential covariance function
#' (\code{ptide_squared_exponential_covariance_function} in this package has
#' the expression
#' \deqn{k(t, t') = \exp\left( -\frac{(t - t')^2}{2(\ell / d)^2} \right),}
#' where \eqn{\ell} is the length scale parameter, and d is a user-specified
#' length dilation.
#' @section Periodic:
#' The periodic covariance function (\code{ptide_periodic_covariance_function}
#' in this package has the expression
#' \deqn{k(t, t')
#'   = \exp\left( -\frac{2\sin^2(\pi|t - t'| / p)}{\ell^2} \right),}
#' where \eqn{\ell} is the length scale parameter, and \eqn{p} is the period of
#' the covariance function.
#' @section Quasi-periodic:
#' The quasi-periodic covariance function
#' (\code{ptide_quasi_periodic_covariance_function} in this package has the
#' expression
#' \deqn{k(t, t') =
#'     \exp\left( -\frac{2\sin^2(\pi|t - t'| / p)}{\ell_p^2} \right)
#'     \exp\left( -\frac{(t - t')^2}{2(\ell_d / d)^2} \right),
#' }
#' where \eqn{\ell_p} is the periodic length scale parameter, \eqn{p} is the
#' period of the covariance function, \eqn{\ell_d} is the decay length scale,
#' and \eqn{d} is the length dilation, as described in the section on the
#' squared exponential covariance function.
#' @section White noise:
#' The white noise covariance function
#' (\code{ptide_white_noise_covariance_function} in this package has the
#' expression
#' \deqn{k(t, t') = I(t = t'),}
#' where \eqn{I(\cdot)} is an indicator function.
#' @param length_dilation Used to multiplicatively scale the hours in the
#' covariance calculation, intended to help the user shift the value of the
#' length scale parameters to approximately one.
#' @param period The period of the periodic portion of a periodic or
#' quasi-periodic covariance function, in hours. Can also be given as a string,
#' in which case the period corresponds to one of the tidal frequencies given
#' in \code{\link{tidal_constituents}}.
#' @param ell_a The shape parameter of the inverse-gamma prior on the length
#' scale parameter, \eqn{\ell}.
#' @param ell_b The scale parameter of the inverse-gamma prior on the length
#' scale parameter, \eqn{\ell}.
#' @param ell_decay_a The shape parameter of the inverse-gamma prior on the
#' decay length scale parameter, ell_decay, of the quasi-periodic covariance
#' function.
#' @param ell_decay_b The rate parameter of the inverse-gamma prior on the decay
#' length scale parameter, ell_decay, of the quasi-periodic covariance function.
#' @param ell_periodic_a The shape parameter of the inverse-gamma prior on the
#' periodic length scale parameter, ell_periodic, of the quasi-periodic
#' covariance function.
#' @param ell_periodic_b The rate parameter of the inverse-gamma prior on the
#' periodic length scale parameter, ell_periodic, of the quasi-periodic
#' covariance function.
#' @param L_shape The parameter of the LKJ prior for the correlation matrix C of
#' the covariance function.
#' @param sigma_squared_a The shape parameter of inverse-gamma prior on sigma
#' squared. If set to 'empirical', this will be based on the empirical variance
#' of the data.
#' @param sigma_squared_b The rate parameter of inverse-gamma prior on sigma
#' squared.
#' @name ptide_covariance_function
NULL

.validate_covariance_function <- function(
  sigma_squared_a,
  sigma_squared_b
) {
  if (!(sigma_squared_a == 'empirical' && sigma_squared_b == 'empirical')) {
    if (!(is.numeric(sigma_squared_a) && is.numeric(sigma_squared_b))) {
      stop('sigma_squared_a and sigma_squared_b must be both numeric or both "empirical"')
    }
  }
}

.covariance_function_parameter_names <- list(
  exponential = 'ell',
  matern32 = 'ell',
  matern52 = 'ell',
  squared_exponential = 'ell',
  periodic = 'ell',
  quasi_periodic = c('ell_periodic', 'ell_decay'),
  white_noise = NULL
)
.covariance_function_configuration_names <- list(
  exponential = 'length_dilation',
  matern32 = 'length_dilation',
  matern52 = 'length_dilation',
  squared_exponential = 'length_dilation',
  periodic = 'period',
  quasi_periodic = c('period', 'length_dilation'),
  white_noise = NULL
)
.covariance_function_types <- names(.covariance_function_parameter_names)
.covariance_function_hyperparameter_names <- list(
  exponential = c(
    'sigma_squared_a',
    'sigma_squared_b',
    'L_shape',
    'ell_a',
    'ell_b'
  ),
  squared_exponential = c(
    'sigma_squared_a',
    'sigma_squared_b',
    'L_shape',
    'ell_a',
    'ell_b'
  ),
  matern32 = c(
    'sigma_squared_a',
    'sigma_squared_b',
    'L_shape',
    'ell_a',
    'ell_b'
  ),
  matern52 = c(
    'sigma_squared_a',
    'sigma_squared_b',
    'L_shape',
    'ell_a',
    'ell_b'
  ),
  periodic = c(
    'period',
    'sigma_squared_a',
    'sigma_squared_b',
    'L_shape',
    'ell_a',
    'ell_b'
  ),
  quasi_periodic = c(
    'period',
    'sigma_squared_a',
    'sigma_squared_b',
    'L_shape',
    'ell_periodic_a',
    'ell_periodic_b',
    'ell_decay_a',
    'ell_decay_b'
  ),
  white_noise = c(
    'sigma_squared_a',
    'sigma_squared_b',
    'L_shape'
  )
)

#' @describeIn ptide_covariance_function Exponential covariance function
#' @export
ptide_exponential_covariance_function <- function(
  sigma_squared_a = 'empirical',
  sigma_squared_b = 'empirical',
  L_shape = 1.1,
  ell_a = 0.1,
  ell_b = 0.1,
  length_dilation = 1
) {
  .validate_covariance_function(sigma_squared_a, sigma_squared_b)
  list(
    type = 'exponential',
    sigma_squared_a = sigma_squared_a,
    sigma_squared_b = sigma_squared_b,
    L_shape = L_shape,
    ell_a = ell_a,
    ell_b = ell_b,
    length_dilation = length_dilation
  )
}

#' @describeIn ptide_covariance_function Matern 3/2 covariance function
#' @export
ptide_matern32_covariance_function <- function(
  sigma_squared_a = 'empirical',
  sigma_squared_b = 'empirical',
  L_shape = 1.1,
  ell_a = 0.1,
  ell_b = 0.1,
  length_dilation = 1
) {
  .validate_covariance_function(sigma_squared_a, sigma_squared_b)
  list(
    type = 'matern32',
    sigma_squared_a = sigma_squared_a,
    sigma_squared_b = sigma_squared_b,
    L_shape = L_shape,
    ell_a = ell_a,
    ell_b = ell_b,
    length_dilation = length_dilation
  )
}

#' @describeIn ptide_covariance_function Matern 5/2 covariance function
#' @export
ptide_matern52_covariance_function <- function(
  sigma_squared_a = 'empirical',
  sigma_squared_b = 'empirical',
  L_shape = 1.1,
  ell_a = 0.1,
  ell_b = 0.1,
  length_dilation = 1
) {
  .validate_covariance_function(sigma_squared_a, sigma_squared_b)
  list(
    type = 'matern52',
    sigma_squared_a = sigma_squared_a,
    sigma_squared_b = sigma_squared_b,
    L_shape = L_shape,
    ell_a = ell_a,
    ell_b = ell_b,
    length_dilation = length_dilation
  )
}

#' @describeIn ptide_covariance_function Periodic covariance function
#' @export
ptide_periodic_covariance_function <- function(
  period = 'M2',
  sigma_squared_a = 'empirical',
  sigma_squared_b = 'empirical',
  L_shape = 1.1,
  ell_a = 0.1,
  ell_b = 0.1,
  length_dilation = 1
) {
  .validate_covariance_function(sigma_squared_a, sigma_squared_b)
  list(
    type = 'periodic',
    period = tidal_period(period),
    sigma_squared_a = sigma_squared_a,
    sigma_squared_b = sigma_squared_b,
    L_shape = L_shape,
    ell_a = ell_a,
    ell_b = ell_b,
    length_dilation = length_dilation
  )
}

#' @describeIn ptide_covariance_function Quasi-periodic covariance function
#' @export
ptide_quasi_periodic_covariance_function <- function(
  period = 'M2',
  sigma_squared_a = 'empirical',
  sigma_squared_b = 'empirical',
  L_shape = 1.1,
  ell_periodic_a = 0.1,
  ell_periodic_b = 0.1,
  ell_decay_a = 0.1,
  ell_decay_b = 0.1,
  length_dilation = 1
) {
  .validate_covariance_function(sigma_squared_a, sigma_squared_b)
  list(
    type = 'quasi_periodic',
    period = tidal_period(period),
    sigma_squared_a = sigma_squared_a,
    sigma_squared_b = sigma_squared_b,
    L_shape = L_shape,
    ell_periodic_a = ell_periodic_a,
    ell_periodic_b = ell_periodic_b,
    ell_decay_a = ell_decay_a,
    ell_decay_b = ell_decay_b,
    length_dilation = length_dilation
  )
}

#' @describeIn ptide_covariance_function White noise covariance function
#' @export
ptide_white_noise_covariance_function <- function(
  sigma_squared_a = 'empirical',
  sigma_squared_b = 'empirical',
  L_shape = 1.1
) {
  .validate_covariance_function(sigma_squared_a, sigma_squared_b)
  list(
    type = 'white_noise',
    sigma_squared_a = sigma_squared_a,
    sigma_squared_b = sigma_squared_b,
    L_shape = L_shape
  )
}
