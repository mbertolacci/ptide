#' Interpolate time series data
#'
#' Interpolate data to make a time series regular and to fill in missing gaps.
#'
#' @param data Input data frame
#' @param fields Field names to interpolate
#' @param times Times to interpolate at
#' @param time_field Name of the time field in \code{data}
#' @param method Interpolation method. Either \code{\link[stats]{approx}} or
#' \code{\link[stats]{spline}}.
#' @param ... Additional arguments passed to the interpolation function.
#'
#' @return A data frame with a time field and each of the fields named in 
#' \code{fields}.
#' @export
interpolate_time_series <- function(
  data,
  fields,
  times,
  time_field = 'time',
  method = c('approx', 'spline'),
  ...
) {
  method <- match.arg(method)

  start_time <- min(data[[time_field]])
  to_minutes <- function(x) as.double(difftime(
    x,
    start_time,
    units = 'mins'
  ))
  input_minutes <- to_minutes(data[[time_field]])
  output_minutes <- to_minutes(times)

  output <- data.frame(
    time = times
  )
  colnames(output)[1] <- time_field
  for (field_name in fields) {
    output[[field_name]] <- match.fun(method)(
      x = input_minutes,
      y = data[[field_name]],
      xout = output_minutes,
      ...
    )$y
  }
  output
}

#' Fill gaps in time series data
#'
#' Fill gaps in time series data by interpolating between known values.
#'
#' @param data Input data frame
#' @param fields Field names to interpolate
#' @param time_field Name of the time field in \code{data}
#' @param delta Time step of the output data. If the input data is regular,
#' this should match the input time step. The special value 'guess' will
#' pick the median time gap; this should not be used with irregular input data.
#' @param max_gap Maximum gap length to fill, in multiples of \code{delta}; if
#' set to \code{Inf}, all gaps will be filled
#' @param ... Additional arguments passed to \code{\link{interpolate_time_series}},
#' which does the gap filling.
#'
#' @return A data frame with a time field and each of the fields named in
#' \code{fields}.
#' @export
fill_time_series_gaps <- function(
  data,
  fields,
  time_field = 'time',
  delta = 'guess',
  max_gap = Inf,
  ...
) {
  if (!is.integer(max_gap) && !is.numeric(max_gap)) {
    stop('max_gap must be a number')
  }
  for (field in fields) {
    if (any(is.na(data[[field]]))) {
      stop(sprintf(
        'Missing values found in "%s"; cannot fill gaps',
        field
      ))
    }
  }

  input_times <- data[[time_field]]
  if (delta == 'guess') {
    delta <- median(diff(input_times))
  }

  output_times <- seq(min(input_times), max(input_times), by = delta)
  output <- interpolate_time_series(
    data,
    fields,
    output_times,
    time_field = time_field,
    ...
  )

  included_rle <- rle(output_times %in% input_times)
  included_rle$values[
    !included_rle$values & included_rle$lengths > max_gap
  ] <- NA
  long_gaps <- is.na(inverse.rle(included_rle))
  for (field in fields) {
    output[[field]][long_gaps] <- NA
  }
  output
}
