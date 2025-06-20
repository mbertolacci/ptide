#' Create a tidal harmonic design matrix
#'
#' Creates a data.frame with harmonics at the chosen frequencies.
#'
#' @param times The times at which to create the harmonics. A vector of date
#' time objects for which difftime is defined.
#' @param hours A numeric vector of times at which to create the harmonics,
#' in hours. Ignored if \code{times} are provided.
#' @param reference_time The reference time (time zero) for the harmonics.
#' Ignored if \code{hours} is provided
#' @param constituent_names The names of the tidal harmonics to include. See
#' \code{\link[ptide]{tidal_constituents}}.
#' @param intercept Whether to include an intercept term as a column of ones.
#' @param trend Whether to include a trend column in hours.
#' @return A data.frame with the chosen frequencies as columns.
#' @examples
#' design_matrix <- tidal_harmonic_design_matrix(hours = 1 : 30)
#' print(design_matrix)
#' @export
tidal_harmonic_design_matrix <- function(
  times,
  hours,
  reference_time = min(times),
  constituent_names = c('M2', 'S2'),
  intercept = FALSE,
  trend = FALSE
) {
  stopifnot(all(constituent_names %in% ptide::tidal_constituents$name))

  constituents <- ptide::tidal_constituents[
    ptide::tidal_constituents$name %in% constituent_names,
    c('name', 'freq')
  ]

  if (!missing(times)) {
    hours <- as.double(difftime(
      times,
      reference_time,
      units = 'hours'
    ))
  }

  do.call(cbind, c(
    if (intercept) data.frame(intercept = rep(1, length(hours))) else NULL,
    if (trend) data.frame(hours = hours) else NULL,
    lapply(1 : nrow(constituents), function(i) {
      constituent <- constituents[i, ]
      output <- cbind(
        sinpi(2 * hours * constituent$freq),
        cospi(2 * hours * constituent$freq)
      )
      colnames(output) <- sprintf(c('%s_sin', '%s_cos'), constituent$name)
      as.data.frame(output)
    })
  ))
}

#' Tidal period
#'
#' Converts a tidal constituent name to its period in hours.
#'
#' @param period A tidal constituent name, or a numeric value in hours.
#' @return The period in hours.
#' @export
tidal_period <- function(period) {
  if (!is.character(period)) return(period)
  1 / (
    ptide::tidal_constituents %>%
      dplyr::filter(name == period) %>%
      dplyr::pull(freq)
  )
}

#' Tidal constituents
#'
#' A data.frame of tidal constituents giving names, frequencies in hours, and
#' Doodson numbers. Copied from \link[oce]{tidedata} in the oce package, who in
#' which in turn derives it from Foreman (1979).
#'
#' @format A data.frame
#' @docType data
'tidal_constituents'

