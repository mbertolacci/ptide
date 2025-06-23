# ptide

## Installation

First, install the remotes package:

```r
install.packages('remotes')
```

Then, install the package:

```r
remotes::install_github('mbertolacci/ptide')
```

## Usage

### Model fitting

Here is a minimal example of how to fit a model to some data.

```r
library(ptide)

# Here we assume that the data is in a data frame called 'observed_df' with
# columns 'easting' and 'northing', and a column 'time' of type POSIXct.

# First you define a ptide_model object, which contains the model specification.
# In this case, we are using a white noise covariance function, and we are
# using the M2, S2, K1, and O1 tidal constituents. This corresponds to a
# classical tidal model, with correlated noise between the easting and northing
# components.
model <- ptide_model(
  cbind(easting, northing) ~ 1,
  tidal_constituents = c('M2', 'S2', 'K1', 'O1'),
  covariance_functions = list(
    ptide_white_noise_covariance_function()
  )
)

# Fit the model; this is pretty fast for the white noise case
fit <- ptide_fit(model, observed_df)

summary(fit)
```

### Predictions

```r
library(ptide)
library(lubridate)

## Assume that the data is in a data frame with columns 'easting' and 'northing',
## and a column 'time' of type POSIXct.

# Construct a data frame of prediction times
predicted_df <- data.frame(
  time = max(observed_df$time) + hours(1 : 24)
)

## Assume that fit is a ptide_fit object

# Make predictions
predictions <- ptide_predict(
  fit,
  predicted_df,
  observed_df
)

# The predictions are in predictions$mean
str(predictions$mean)
```
