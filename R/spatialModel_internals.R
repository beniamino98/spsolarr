#' Compute monthly moments for mixture with 16 components
#' @param nb_models a list with 4 neighborhoods models
#' @param nmonths numeric vector, reference months. Default is `1:12`, i.e. all the months.
#' @param weights numeric, vector of 4 weights.
#' @noRd
#' @keywords spatialModel internal
spatialModel_combinations = function(nb_models, nmonths = 1:12, weights = rep(1/4, 4), nobs.min = 5){
  # Weights in vector notation ensuring normalization
  w <- matrix(weights/sum(weights), ncol = 1)
  # Possible combinations in a grid with 4 reference locations
  possible_combinations <- data.frame(
    State = 1:16,
    Month = rep(1, 16),
    Location_A = c(1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0),
    Location_B = c(1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0),
    Location_C = c(1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0),
    Location_D = c(1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0)
  )
  # Add observations weights
  data_A <- dplyr::bind_cols(nb_models[[1]]$data)
  data_B <- dplyr::bind_cols(nb_models[[2]]$data)
  data_C <- dplyr::bind_cols(nb_models[[3]]$data)
  data_D <- dplyr::bind_cols(nb_models[[4]]$data)
  # Filter only for train data with non-zero weight
  df_A <- dplyr::filter(data_A, isTrain & weights != 0)
  df_B <- dplyr::filter(data_B, isTrain & weights != 0)
  df_C <- dplyr::filter(data_C, isTrain & weights != 0)
  df_D <- dplyr::filter(data_D, isTrain & weights != 0)
  # Select only relevant variables
  df_A <- dplyr::select(df_A, date, B_A = "B", u_A = "u_tilde")
  df_B <- dplyr::select(df_B, date, B_B = "B", u_B = "u_tilde")
  df_C <- dplyr::select(df_C, date, B_C = "B", u_C = "u_tilde")
  df_D <- dplyr::select(df_D, date, B_D = "B", u_D = "u_tilde")
  # Create a unique dataset
  data <- dplyr::left_join(df_A, df_B, by = "date") %>%
    dplyr::left_join(df_C, by = "date") %>%
    dplyr::left_join(df_D, by = "date") %>%
    dplyr::mutate(Month = lubridate::month(date))

  # Main loop
  m <- 1
  combinations <- list()
  for(m in nmonths){
    # Initialize the dataset
    combinations[[m]] <- dplyr::as_tibble(possible_combinations)
    combinations[[m]]$probs <- 0
    combinations[[m]]$mean <- 0
    combinations[[m]]$sd <- 0
    combinations[[m]]$n <- 0
    combinations[[m]]$cv <- as.list(1:16)
    # Filter dataset for reference month
    df <- dplyr::filter(data, Month == m)
    s <- 1
    for(s in 1:16){
      # Extract the states
      states <- unlist(possible_combinations[s,][-c(1:2)])
      # Compute joint probabilities
      combinations[[m]]$probs[s] <- mean(df$B_A == states[1] & df$B_B == states[2] & df$B_C == states[3] & df$B_D == states[4], na.rm = TRUE)
      # Extract the realized joint series
      df_u <- dplyr::filter(df, B_A == states[1] & B_B == states[2] & B_C == states[3] & B_D == states[4])
      df_u <- dplyr::select(df_u, dplyr::contains("u"))
      # Number of observations for the estimates
      combinations[[m]]$n[s] <- nrow(df_u)
      # If the number of cases is too low set probability to zero
      if (nrow(df_u) <= nobs.min) {
        combinations[[m]]$probs[s] <- 0
        combinations[[m]]$sd[s] <- 1
        next
      }
      # Variance-covariance matrix
      combinations[[m]]$cv[[s]] <- cov(df_u)
      # Expected value for the i-th state
      combinations[[m]]$mean[s] <- sum(dplyr::summarise_all(df_u, mean)*w)
      combinations[[m]]$sd[s] <- sqrt(t(w) %*% combinations[[m]]$cv[[s]] %*% w)[1]
    }
    # Normalize the probabilities
    combinations[[m]]$probs <- combinations[[m]]$probs/sum(combinations[[m]]$probs)
    combinations[[m]]$Month <- m
    colnames(combinations[[m]]) <- c("State", "Month", names(nb_models), "probs", "mean", "sd","n", "cv", "mu2")
  }
  dplyr::bind_rows(combinations)
}

#' Perform the bilinear interpolation for a target variable.
#' @param nb_models a list with 4 neighborhoods models
#' @param target character, name of the target variable to interpolate.
#' @param n number of neighborhoods to use for interpolation.
#' @param weights numeric, vector of weights of with length equal to n.
#' @param day_date date for interpolation, if missing all the available dates will be used.
#' @noRd
#' @keywords spatialModel internal
spatialModel_interpolator = function(nb_models, target = "GHI", n = 4, weights = rep(1/n, n), day_date){
  # Initialize the output dataset
  interp_data <- dplyr::tibble(x3 = lubridate::as_date(NA), x4 = NA, x5 = FALSE)
  colnames(interp_data) <- c("date", target, "interpolated")
  # Neighborhoods models
  interp_data <- as.list(interp_data)
  if (!missing(day_date)) {
    data <- purrr::map(nb_models, ~dplyr::filter(.x$data, date %in% as.Date(day_date)))
    data <- purrr::map(data, ~.x[, c("date", target)])
  } else {
    data <- purrr::map(nb_models, ~.x$data[, c("date", target)])
  }
  df <- data[[1]]
  for(i in 2:n){
    df <- dplyr::left_join(df, data[[i]], by = "date")
  }
  interp_data$date <- df$date
  interp_data[[target]] <- 0
  for(i in 1:n){
    interp_data[[target]] <- interp_data[[target]] + df[,i+1][[1]]*weights[i]
  }
  interp_data <- dplyr::bind_cols(interp_data)
  interp_data$interpolated <- TRUE
  return(interp_data)
}
