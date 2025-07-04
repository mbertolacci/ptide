% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{gamma_quantile_prior}
\alias{gamma_quantile_prior}
\alias{inverse_gamma_quantile_prior}
\title{Construct Gamma and Inverse Gamma Distributions with Chosen Quantiles}
\usage{
gamma_quantile_prior(
  q_lower,
  q_upper,
  p_lower = 0.05,
  p_upper = 0.95,
  interval = c(1e-04, 1000),
  tolerance = sqrt(.Machine$double.eps)
)

inverse_gamma_quantile_prior(q_lower, q_upper, ...)
}
\arguments{
\item{q_lower}{Lower quantile value}

\item{q_upper}{Upper quantile value}

\item{p_lower}{Probability corresponding to the lower quantile.}

\item{p_upper}{Probability corresponding to the upper quantile.}

\item{interval}{Interval in which to search for the shape parameter.}

\item{tolerance}{Tolerance for the algorithm.}

\item{...}{Passed from \code{inverse_gamma_quantile_prior} to
\code{gamma_quantile_prior}.}
}
\value{
A list containing the appropriate shape and rate parameters for the
distribution.
}
\description{
These functions construct a gamma/inverse gamma distribution by finding the
shape and rate parameters that match the given quantiles and probabilities.
It can be useful for constructing informative gamma/inverse gamma prior
distributions.
}
\section{Functions}{
\itemize{
\item \code{inverse_gamma_quantile_prior()}: Equivalent function for the inverse gamma
distribution.

}}
\examples{
# Create a gamma distribution that has 5th and 95th percentiles at 2 and 10
# respectively
gamma_params <- gamma_quantile_prior(q_lower = 2, q_upper = 10)
print(gamma_params)
}
\seealso{
\code{\link[stats]{qgamma}}
}
