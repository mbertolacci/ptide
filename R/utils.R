#' Construct Gamma and Inverse Gamma Distributions with Chosen Quantiles
#'
#' These functions construct a gamma/inverse gamma distribution by finding the
#' shape and rate parameters that match the given quantiles and probabilities.
#' It can be useful for constructing informative gamma/inverse gamma prior
#' distributions.
#'
#' @param q_lower Lower quantile value
#' @param q_upper Upper quantile value
#' @param p_lower Probability corresponding to the lower quantile.
#' @param p_upper Probability corresponding to the upper quantile.
#' @param interval Interval in which to search for the shape parameter.
#' @param tolerance Tolerance for the algorithm.
#' @param ... Passed from \code{inverse_gamma_quantile_prior} to
#' \code{gamma_quantile_prior}.
#' @return A list containing the appropriate shape and rate parameters for the
#' distribution.
#' @examples
#' # Create a gamma distribution that has 5th and 95th percentiles at 2 and 10
#' # respectively
#' gamma_params <- gamma_quantile_prior(q_lower = 2, q_upper = 10)
#' print(gamma_params)
#' @seealso
#' \code{\link[stats]{qgamma}}
#' @export
gamma_quantile_prior <- function(
  q_lower, q_upper,
  p_lower = 0.05, p_upper = 0.95,
  interval = c(0.0001, 1000),
  tolerance = sqrt(.Machine$double.eps)
) {
  stopifnot(q_lower > 0)
  stopifnot(q_lower < q_upper)
  stopifnot(p_lower < p_upper)

  # Find the shape parameter that gives the appropriate ratio between the
  # quantiles
  ratio <- q_lower / q_upper
  shape <- uniroot(function(shape) {
    theoretical <- qgamma(c(p_lower, p_upper), shape = shape, rate = 1)
    theoretical[1] / theoretical[2] - ratio
  }, interval, tol = tolerance)$root

  # Find the rate parameter that gives the correct quantiles
  rate <- qgamma(p_upper, shape = shape, rate = 1) / q_upper

  list(shape = shape, rate = rate)
}

#' @describeIn gamma_quantile_prior Equivalent function for the inverse gamma
#' distribution.
#' @export
inverse_gamma_quantile_prior <- function(
  q_lower, q_upper, ...
) {
  gamma_quantile_prior(1 / q_upper, 1 / q_lower, ...)
}


# As per http://stackoverflow.com/questions/14500707/select-along-one-of-n-dimensions-in-array
.index_array <- function(x, dim, i, drop = FALSE) {
  # Create list representing arguments supplied to [
  # bquote() creates an object corresponding to a missing argument
  indices <- rep(list(bquote()), length(dim(x)))
  indices[[dim]] <- i

  # Generate the call to [
  call <- as.call(c(
    list(as.name('['), quote(x)),
    indices,
    list(drop = drop)
  ))
  # Finally, evaluate it
  eval(call)
}

'.index_array<-' <- function(x, dim, i, value) {
  # Create list representing arguments supplied to [
  # bquote() creates an object corresponding to a missing argument
  indices <- rep(list(bquote()), length(dim(x)))
  indices[[dim]] <- i

  # Generate the call to [
  call <- as.call(c(
    list(as.name('[<-'), quote(x)),
    indices,
    list(quote(value))
  ))
  # Finally, evaluate it
  eval(call)
}
