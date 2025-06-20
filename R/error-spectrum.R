#' Error spectrum
#'
#' Calculate (and optionally plot) the spectrum of the error process for a ptide
#' fit. The spectrum is calculated through numerical integration based on the
#' discrete Fourier transform. The estimated spectral density matrix is returned
#' for a range of frequencies determined by \code{frequency_range},
#' \code{delta}, and \code{lag_max}.
#'
#' The accuracy of the returned spectrum depends on both the \code{delta} and
#' \code{lag_max} arguments. These should be adjusted (made smaller for
#' \code{delta} and larger for \code{lag_max}) until the spectrum converges.
#'
#' @param fit A ptide fit object made using \code{\link{ptide_fit}}.
#' @param frequency_range A two-element vector giving the range of frequencies
#' to return. These are measured in cycles per hour.
#' @param delta The time step used in the numerical integration. This should be
#' made as small as possible subject to computing time.
#' @param lag_max The maximum lag to use in the numerical integration. This
#' should be large enough that the autocovariance is essentially zero at this
#' lag.
#' @param plot Whether to plot the spectrum.
#' @param ... Additional arguments passed to \code{\link{plot}}.
#'
#' @return An n by d by d array. The first dimension corresponds to the
#' frequencies. The second and third dimensions correspond to the rows and
#' columns of the spectral matrix. The frequencies are stored as an attribute
#' named \code{frequencies}.
#' @export
ptide_error_spectrum <- function(
  fit,
  frequency_range = c(0, 1),
  delta = 0.1,
  lag_max = 500,
  plot = FALSE,
  ...
) {
  max_frequency <- 0.5 / delta
  if (frequency_range[2] > max_frequency) {
    stop('frequency_range[2] must be <= 0.5 / delta')
  }
  lags <- seq(0, lag_max, by = delta)
  n_lags <- length(lags)
  acf_lags <- ptide_cross_covariance(lags, 0, fit)

  d <- dim(acf_lags)[3]
  acf_lags_matrix <- matrix(acf_lags, nrow = n_lags)
  fft_lags_matrix <- mvfft(acf_lags_matrix)

  fft_frequencies <- (0 : (n_lags - 1)) / (delta * n_lags)
  output_indices <- (
    fft_frequencies >= frequency_range[1]
    & fft_frequencies <= frequency_range[2]
  )
  spectrum_matrix <- 2 * delta * Re(fft_lags_matrix)[output_indices, , drop = FALSE]

  output <- array(
    spectrum_matrix,
    dim = dim(acf_lags[output_indices, 1, , ])
  )
  attr(output, 'frequencies') <- fft_frequencies[output_indices]

  if (plot) {
    if (d == 1) {
      plot(
        attr(output, 'frequencies'),
        output,
        type = 'l',
        xlab = 'Frequency [1 / hour]',
        ylab = 'Spectrum',
        ...
      )

    } else {
      withr::with_par(list(mfrow = c(d, d)), {
        for (i in seq_len(d)) {
          for (j in seq_len(d)) {
            plot(
              attr(output, 'frequencies'),
              output[, i, j],
              type = 'l',
              xlab = 'Frequency [1 / hour]',
              ylab = 'Spectrum',
              main = sprintf('(%d, %d)', i, j),
              ...
            )
          }
        }
      })
    }
    invisible(output)
  } else {
    output
  }
}
