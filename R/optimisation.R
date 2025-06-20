#' Find Optimal Parameters for a Stan Model Using nlminb
#'
#' This function optimizes the parameters of a given Stan model using
#' \code{\link{nlminb}}. This function internally uses the
#' \code{\link[rstan]{log_prob}} and \code{\link[rstan]{grad_log_prob}}
#' functions from the \pkg{rstan} package.
#'
#' @param object A Stan model object.
#' @param data A list containing the data to be used in the Stan model.
#' @param trace How frequently to print the current iteration (every `trace`
#' iterations). Default is 2.
#' @param control A list of control parameters passed to
#' \code{\link{nlminb}}.
#' @param initial_value_attempts The number of attempts to find a valid initial
#'   value for the optimisation. Default is 10.
#' @param theta_initial The initial value for the parameters. If not supplied,
#' a random value will be used.
#' @param ... Additional arguments passed to \code{\link{nlminb}}.
#'
#' @return A list containing the following:
#' \itemize{
#'   \item \code{value}: The optimized value of the objective function.
#'   \item \code{par}: The optimized parameters in the original parameter space.
#'   \item \code{iterations}: The number of iterations taken by the algorithm.
#'   \item \code{evaluations}: The number of function evaluations.
#'   \item \code{convergence}: The convergence code returned by \code{\link{nlminb}}.
#'   \item \code{message}: The message returned by \code{\link{nlminb}}.
#' }
#'
#' @examples
#' \dontrun{
#'   # Define a Stan model (assuming `model_code` contains your Stan model)
#'   stan_model <- rstan::stan_model(model_code = model_code)
#'   # Define data (assuming `data_list` contains your data)
#'   fit <- stan_nlminb(stan_model, data = data_list)
#' }
#'
#' @seealso
#' \code{\link{nlminb}}
#' @export
stan_nlminb <- function(
  object,
  data,
  trace = 2,
  control = list(
    iter.max = 2000,
    eval.max = 10000,
    trace = trace
  ),
  initial_value_attempts = 10,
  theta_initial,
  ...
) {
  fit0 <- suppressMessages(rstan::sampling(object, data, chain = 0))
  n_params <- rstan::get_num_upars(fit0)

  if (is.null(control$trace)) {
    control$trace <- trace
  }

  if (missing(theta_initial)) {
    for (attempt in seq_len(initial_value_attempts)) {
      theta_initial <- runif(n_params, -4, 4)
      value <- try(rstan::log_prob(fit0, theta_initial))
      if (!inherits(value, 'try-error') && is.finite(value)) {
        break
      }
    }
    if (attempt == initial_value_attempts) {
      stop('Failed to find a valid initial value')
    }
  }

  fit <- nlminb(
    theta_initial,
    function(theta) {
      tryCatch({
        -rstan::log_prob(fit0, theta, adjust_transform = FALSE)
      }, error = function(e) {
        Inf
      })
    },
    function(theta) {
      -rstan::grad_log_prob(fit0, theta, adjust_transform = FALSE)
    },
    control = control,
    ...
  )

  list(
    value = -fit$objective,
    par = rstan::constrain_pars(fit0, fit$par),
    iterations = fit$iterations,
    evaluations = fit$evaluations,
    convergence = fit$convergence,
    message = fit$message
  )
}

#' Find Optimal Parameters for a Stan Model Using Multiple Attempts
#'
#' This function optimizes the parameters of a given Stan model using multiple
#' attempts by calling either \code{\link[rstan]{optimizing}} or
#' \code{\link{stan_nlminb}}. The function returns the best result of
#' `.best_of` attempts. It will attempt up to `.max_attempts` times to find a
#' solution, which makes it somewhat robust to errors in the optimisation
#' process. Logging in this function uses the \pkg{logger} package, and can be
#' enabled by setting the log level to trace through
#' `logger::log_threshold(logger::TRACE)`.
#'
#' @param ... Arguments passed to the chosen method.
#' @param .max_attempts The maximum number of attempts. Default is 20.
#' @param .best_of The number of attempts to keep. Default is 10.
#' @param .method The optimisation function tocall. Default is
#' \code{\link[rstan]{optimizing}}; can also be \code{stan_trust_optim}. Any
#' function with return type compliant with `rstan::optimizing` can be used.
#'
#' @return A list containing the result of the best attempt and all attempts.
#'
#' @examples
#' \dontrun{
#'   # Define a Stan model (assuming `model_code` contains your Stan model)
#'   stan_model <- rstan::stan_model(model_code = model_code)
#'   # Define data (assuming `data_list` contains your data)
#'   fit <- stan_best_of_optim(rstan::optimizing, stan_model, data = data_list)
#' }
#'
#' @export
stan_optimizing_best_of <- function(
  ...,
  .max_attempts = 20,
  .best_of = 10,
  .method = rstan::optimizing
) {
  .method <- match.fun(.method)

  good_attempts <- 0
  best_result <- list(value = -Inf)
  attempts <- list()
  for (attempt in seq_len(.max_attempts)) {
    logger::log_trace('Running attempt {attempt} of {.max_attempts} | completed {good_attempts} of {.best_of} | best {best_result$value}')
    result <- tryCatch({
      .method(...)
    }, error = function(e) {
      print(e)
      NULL
    })
    if (!is.null(result)) {
      attempts <- c(attempts, list(result))
    }
    if (is.null(result)) next

    good_attempts <- good_attempts + 1
    if (result$value > best_result$value) {
      best_result <- result
    }
    if (good_attempts == .best_of) break
  }

  if (attempt == .max_attempts) stop('Max attempts exceeded')

  best_result$attempts <- attempts
  best_result
}
