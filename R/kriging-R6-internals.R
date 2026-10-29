#' Cross-validate multiple kriging specifications
#'
#' Compares several simple, ordinary, and universal kriging specifications by
#' leave-one-out cross-validation (LOOCV). The function evaluates all requested
#' combinations of correlation models and, for universal kriging, trend formulas,
#' then returns the models ranked by cross-validation negative log-likelihood.
#'
#' @param x Numeric vector of observed values to be kriged.
#' @param lat Numeric vector of latitudes associated with \code{x}.
#' @param lon Numeric vector of longitudes associated with \code{x}.
#' @param cor_formulas Optional character vector of trend formulas to test in
#'   universal kriging. If missing, a default set of polynomial trend formulas
#'   in latitude and longitude is used.
#' @param cor_models Optional character vector of spatial correlation models to
#'   test. If missing, defaults to \code{c("Exp", "Sph", "Mat", "Gau")}.
#' @param n.neighb Integer. Number of nearest neighbours used for prediction in
#'   all cross-validation routines. Default is \code{10}.
#'
#' @details
#' The function evaluates:
#' \itemize{
#'   \item all specified simple kriging models;
#'   \item all specified ordinary kriging models;
#'   \item all specified universal kriging models for each combination of
#'   \code{cor_models} and \code{cor_formulas}.
#' }
#'
#' Results from all candidate specifications are stacked into a single tibble
#' and sorted by \code{CVNLL} in ascending order.
#'
#' @return A tibble with one row per candidate kriging specification and columns:
#' \describe{
#'   \item{\code{kriging}}{Type of kriging model: \code{"sk"}, \code{"ok"}, or \code{"uk"}.}
#'   \item{\code{model}}{Correlation model used.}
#'   \item{\code{formula}}{Trend formula used, or \code{NA} when not applicable.}
#'   \item{\code{CVNLL}}{Mean cross-validation negative log-likelihood.}
#'   \item{\code{CVRMSE}}{Cross-validation root mean squared error.}
#'   \item{\code{theta}}{List-column containing observed values.}
#'   \item{\code{m_theta}}{List-column containing LOOCV predictive means.}
#'   \item{\code{v_theta}}{List-column containing LOOCV predictive variances.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' lat <- runif(30, 44, 45)
#' lon <- runif(30, 11, 12)
#' x   <- 1 + 0.5 * lat - 0.3 * lon + rnorm(30, sd = 0.2)
#'
#' out <- cross_validate_kriging_model(
#'   x = x,
#'   lat = lat,
#'   lon = lon,
#'   cor_models = c("Exp", "Gau"),
#'   cor_formulas = c("x ~ 1", "x ~ 1 + lat + lon"),
#'   n.neighb = 10
#' )
#'
#' out
#' }
#'
#' @export
cross_validate_kriging_model <- function(x, lat, lon, cor_formulas, cor_models, n.neighb = 10){
  # Formulas to test 
  if (missing(cor_formulas)){
    cor_formulas <- c("x ~ 1", 
                      "x ~ 1 + lat + lon", 
                      "x ~ 1 + lat + lon + I(lat*lon)", 
                      "x ~ 1 + lat + lon + I(lat*lon) + I(lat^2) + I(lon^2)", 
                      "x ~ 1 + lat + lon + I(lat^2) + I(lon^2)", 
                      "x ~ 1 + lat + lon + I(lat*lon) + I(lat^2) + I(lon^2) + I(lat^2*lon^2)")
  }
  # Correlation models to test
  if (missing(cor_models)){
    cor_models <- c("Exp", "Sph", "Mat", "Gau")
  }
  
  # Simple kriging
  k <- 1
  sk_models <- list()
  for(cor_model in cor_models){
    print(paste0("simpleKriging: ", cor_model))
    sk_models[[k]] <- cross_validate_simple_kriging(x, lat, lon, cor_model, n.neighb)
    k <- k +1 
  }
  sk_models <- dplyr::bind_cols(kriging = "sk", dplyr::bind_rows(sk_models))
  
  # Ordinary kriging
  k <- 1
  ok_models <- list()
  for(cor_model in cor_models){
    print(paste0("ordinaryKriging: ", cor_model))
    ok_models[[k]] <- cross_validate_ordinary_kriging(x, lat, lon, cor_model, n.neighb)
    k <- k +1 
  }
  ok_models <- dplyr::bind_cols(kriging = "ok", dplyr::bind_rows(ok_models))
  
  # Universal kriging
  k <- 1
  uk_models <- list()
  for(cor_model in cor_models){
    for(cor_formula in cor_formulas){
      print(paste0("universalKriging: ", cor_model, " with formula: ", cor_formula))
      uk_models[[k]] <- cross_validate_universal_kriging(x, lat, lon, cor_formula, cor_model, n.neighb)
      k <- k + 1
    }
  }
  uk_models <- dplyr::bind_cols(kriging = "uk", dplyr::bind_rows(uk_models))
  # All models 
  all_models <- dplyr::bind_rows(sk_models, ok_models, uk_models) 
  dplyr::arrange(all_models, CVNLL)
}
#' Leave-one-out cross-validation for universal kriging
#'
#' Performs leave-one-out cross-validation (LOOCV) for a universal kriging
#' model. At each iteration, one observation is removed, the kriging model is
#' refitted on the remaining locations, and the omitted value is predicted
#' together with its kriging variance.
#'
#' The function returns two cross-validation scores:
#' \itemize{
#'   \item \code{CVNLL}: mean Gaussian negative log-likelihood based on the
#'   predictive mean and variance;
#'   \item \code{CVRMSE}: root mean squared prediction error.
#' }
#'
#' @param x Numeric vector of observed values to be kriged.
#' @param lat Numeric vector of latitudes associated with \code{x}.
#' @param lon Numeric vector of longitudes associated with \code{x}.
#' @param cor_formula Character string specifying the trend formula used in the
#'   universal kriging model. Default is \code{"x ~ 1"}.
#' @param cor_model Character string specifying the spatial correlation model
#'   passed to \code{universalKriging$fit()}. Default is \code{"Exp"}.
#' @param n.neighb Integer. Number of nearest neighbours used for prediction.
#'   Default is \code{10}.
#' @param quiet Logical. If \code{FALSE}, prints the current LOOCV iteration.
#'   Default is \code{FALSE}.
#'
#' @details
#' For each location \eqn{\ell = 1, \dots, L}, the function:
#' \enumerate{
#'   \item removes observation \eqn{\ell};
#'   \item fits a universal kriging model on the remaining \eqn{L-1} locations;
#'   \item predicts the omitted value at \eqn{(\mathrm{lat}_\ell, \mathrm{lon}_\ell)};
#'   \item stores the predictive mean and variance.
#' }
#'
#' The cross-validation negative log-likelihood is computed as
#' \deqn{
#' \mathrm{CVNLL} =
#' \frac{1}{L}\sum_{\ell=1}^L
#' \left[
#' \frac{1}{2}\log(2\pi v_\ell)
#' +
#' \frac{(\theta_\ell - \hat{\theta}_\ell)^2}{2v_\ell}
#' \right],
#' }
#' where \eqn{\theta_\ell} is the true value, \eqn{\hat{\theta}_\ell} the
#' leave-one-out prediction, and \eqn{v_\ell} the corresponding kriging variance.
#'
#' @return A tibble with one row and the following columns:
#' \describe{
#'   \item{\code{model}}{Correlation model used.}
#'   \item{\code{formula}}{Trend formula used in universal kriging.}
#'   \item{\code{CVNLL}}{Mean cross-validation negative log-likelihood.}
#'   \item{\code{CVRMSE}}{Cross-validation root mean squared error.}
#'   \item{\code{theta}}{List-column containing the observed values.}
#'   \item{\code{m_theta}}{List-column containing the LOOCV predictive means.}
#'   \item{\code{v_theta}}{List-column containing the LOOCV predictive variances.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' lat <- runif(20, 44, 45)
#' lon <- runif(20, 11, 12)
#' x   <- 2 + 0.3 * lat - 0.2 * lon + rnorm(20, sd = 0.1)
#'
#' cross_validate_universal_kriging(
#'   x = x,
#'   lat = lat,
#'   lon = lon,
#'   cor_formula = "x ~ 1 + lat + lon",
#'   cor_model = "Exp",
#'   n.neighb = 8
#' )
#' }
#'
#' @export
cross_validate_universal_kriging <- function(x, lat, lon, cor_formula = "x ~ 1", cor_model = "Exp", n.neighb = 10, quiet = FALSE){
  # Full dataset with all the locations
  full_data <- data.frame(x = x, lat = lat, lon = lon)
  L <- nrow(full_data)
  
  # True parameter 
  theta   <- vector("numeric", length = L)
  m_theta <- vector("numeric", length = L)
  v_theta <- vector("numeric", length = L)
  for(l in 1:L) {
    if(!quiet) message(l, "/", L, "\r", appendLF = FALSE)
    # Observation l 
    data_l <- full_data[l, ]
    lat0 <- data_l$lat
    lon0 <- data_l$lon
    theta[l] <- data_l$x[[1]]
    # Remove observation l from full data
    reduced_data <- full_data[-l, ]
    # Fit the model without the l-th parameter 
    uk <- universalKriging$new(reduced_data$x, reduced_data$lat, reduced_data$lon)
    uk$fit(formula = cor_formula,  model = cor_model)
    # Predict mean and variance of the l-th parameter
    theta_hat <- purrr::quietly(uk$predict)(lat = lat0, lon = lon0, n = n.neighb)$result
    m_theta[l] <- theta_hat$theta
    v_theta[l] <- theta_hat$var_theta
  }
  
  # CV score -> best minimize CVNLL 
  CVNLL <- sum(0.5 * log(2 * base::pi * v_theta) + (theta - m_theta)^2 / (2 * v_theta)) / L
  # CV-RMSE (simpler)
  CVRMSE <- sqrt(sum((theta - m_theta)^2) / L)
  
  # Results 
  dplyr::tibble(
    model = cor_model,
    formula = cor_formula,
    CVNLL = CVNLL, 
    CVRMSE = CVRMSE, 
    theta = list(theta), 
    m_theta = list(m_theta),
    v_theta = list(v_theta)
  )
}

#' Leave-one-out cross-validation for simple kriging
#'
#' Performs leave-one-out cross-validation (LOOCV) for a simple kriging model.
#' At each iteration, one observation is removed, the model is refitted on the
#' remaining locations using the sample mean of the reduced dataset as known
#' mean, and the omitted value is predicted together with its kriging variance.
#'
#' The function returns two cross-validation scores:
#' \itemize{
#'   \item \code{CVNLL}: mean Gaussian negative log-likelihood based on the
#'   predictive mean and variance;
#'   \item \code{CVRMSE}: root mean squared prediction error.
#' }
#'
#' @param x Numeric vector of observed values to be kriged.
#' @param lat Numeric vector of latitudes associated with \code{x}.
#' @param lon Numeric vector of longitudes associated with \code{x}.
#' @param cor_model Character string specifying the spatial correlation model
#'   passed to \code{simpleKriging$fit()}. Default is \code{"Exp"}.
#' @param n.neighb Integer. Number of nearest neighbours used for prediction.
#'   Default is \code{10}.
#' @param quiet Logical. If \code{FALSE}, prints the current LOOCV iteration.
#'   Default is \code{FALSE}.
#'
#' @details
#' For each location \eqn{\ell = 1, \dots, L}, the function:
#' \enumerate{
#'   \item removes observation \eqn{\ell};
#'   \item computes the sample mean of the reduced dataset;
#'   \item fits a simple kriging model on the remaining \eqn{L-1} locations;
#'   \item predicts the omitted value at \eqn{(\mathrm{lat}_\ell, \mathrm{lon}_\ell)};
#'   \item stores the predictive mean and variance.
#' }
#'
#' The cross-validation negative log-likelihood is computed assuming Gaussian
#' predictive errors with kriging variance.
#'
#' @return A tibble with one row and the following columns:
#' \describe{
#'   \item{\code{model}}{Correlation model used.}
#'   \item{\code{formula}}{\code{NA}, since no trend formula is used.}
#'   \item{\code{CVNLL}}{Mean cross-validation negative log-likelihood.}
#'   \item{\code{CVRMSE}}{Cross-validation root mean squared error.}
#'   \item{\code{theta}}{List-column containing the observed values.}
#'   \item{\code{m_theta}}{List-column containing the LOOCV predictive means.}
#'   \item{\code{v_theta}}{List-column containing the LOOCV predictive variances.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' lat <- runif(20, 44, 45)
#' lon <- runif(20, 11, 12)
#' x   <- rnorm(20)
#'
#' cross_validate_simple_kriging(
#'   x = x,
#'   lat = lat,
#'   lon = lon,
#'   cor_model = "Gau",
#'   n.neighb = 8
#' )
#' }
#'
#' @export
cross_validate_simple_kriging <- function(x, lat, lon, cor_model = "Exp", n.neighb = 10, quiet = FALSE){
  # Full dataset with all the locations
  full_data <- data.frame(x = x, lat = lat, lon = lon)
  L <- nrow(full_data)
  
  # True parameter 
  theta   <- vector("numeric", length = L)
  m_theta <- vector("numeric", length = L)
  v_theta <- vector("numeric", length = L)
  for(l in 1:L) {
    if(!quiet) message(l, "/", L, "\r", appendLF = FALSE)
    # Observation l 
    data_l <- full_data[l, ]
    lat0 <- data_l$lat
    lon0 <- data_l$lon
    theta[l] <- data_l$x[[1]]
    # Remove observation l from full data
    reduced_data <- full_data[-l, ]
    # Fit the model without the l-th parameter 
    sk <- simpleKriging$new(reduced_data$x, reduced_data$lat, reduced_data$lon)
    sk$fit(mu_theta = mean(reduced_data$x), model = cor_model)
    # Predict mean and variance of the l-th parameter
    theta_hat <- purrr::quietly(sk$predict)(lat = lat0, lon = lon0, n = n.neighb)$result
    m_theta[l] <- theta_hat$theta
    v_theta[l] <- theta_hat$var_theta
  }
  
  # CV score -> best minimize CVNLL 
  CVNLL <- sum(0.5 * log(2 * base::pi * v_theta) + (theta - m_theta)^2 / (2 * v_theta)) / L
  # CV-RMSE (simpler)
  CVRMSE <- sqrt(sum((theta - m_theta)^2) / L)
  
  # Results 
  dplyr::tibble(
    model = cor_model,
    formula = NA,
    CVNLL = CVNLL, 
    CVRMSE = CVRMSE, 
    theta = list(theta), 
    m_theta = list(m_theta),
    v_theta = list(v_theta)
  )
}

#' Leave-one-out cross-validation for ordinary kriging
#'
#' Performs leave-one-out cross-validation (LOOCV) for an ordinary kriging
#' model. At each iteration, one observation is removed, the model is refitted
#' on the remaining locations, and the omitted value is predicted together with
#' its kriging variance.
#'
#' The function returns two cross-validation scores:
#' \itemize{
#'   \item \code{CVNLL}: mean Gaussian negative log-likelihood based on the
#'   predictive mean and variance;
#'   \item \code{CVRMSE}: root mean squared prediction error.
#' }
#'
#' @param x Numeric vector of observed values to be kriged.
#' @param lat Numeric vector of latitudes associated with \code{x}.
#' @param lon Numeric vector of longitudes associated with \code{x}.
#' @param cor_model Character string specifying the spatial correlation model
#'   passed to \code{ordinaryKriging$fit()}. Default is \code{"Exp"}.
#' @param n.neighb Integer. Number of nearest neighbours used for prediction.
#'   Default is \code{10}.
#' @param quiet Logical. If \code{FALSE}, prints the current LOOCV iteration.
#'   Default is \code{FALSE}.
#'
#' @details
#' For each location \eqn{\ell = 1, \dots, L}, the function:
#' \enumerate{
#'   \item removes observation \eqn{\ell};
#'   \item fits an ordinary kriging model on the remaining \eqn{L-1} locations;
#'   \item predicts the omitted value at \eqn{(\mathrm{lat}_\ell, \mathrm{lon}_\ell)};
#'   \item stores the predictive mean and variance.
#' }
#'
#' The cross-validation negative log-likelihood is computed assuming Gaussian
#' predictive errors with kriging variance.
#'
#' @return A tibble with one row and the following columns:
#' \describe{
#'   \item{\code{model}}{Correlation model used.}
#'   \item{\code{formula}}{\code{NA}, since no explicit trend formula is used.}
#'   \item{\code{CVNLL}}{Mean cross-validation negative log-likelihood.}
#'   \item{\code{CVRMSE}}{Cross-validation root mean squared error.}
#'   \item{\code{theta}}{List-column containing the observed values.}
#'   \item{\code{m_theta}}{List-column containing the LOOCV predictive means.}
#'   \item{\code{v_theta}}{List-column containing the LOOCV predictive variances.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' lat <- runif(20, 44, 45)
#' lon <- runif(20, 11, 12)
#' x   <- rnorm(20)
#'
#' cross_validate_ordinary_kriging(
#'   x = x,
#'   lat = lat,
#'   lon = lon,
#'   cor_model = "Mat",
#'   n.neighb = 8
#' )
#' }
#'
#' @export
cross_validate_ordinary_kriging <- function(x, lat, lon, cor_model = "Exp", n.neighb = 10, quiet = FALSE){
  # Full dataset with all the locations
  full_data <- data.frame(x = x, lat = lat, lon = lon)
  L <- nrow(full_data)
  
  # True parameter 
  theta   <- vector("numeric", length = L)
  m_theta <- vector("numeric", length = L)
  v_theta <- vector("numeric", length = L)
  for(l in 1:L) {
    if(!quiet) message(l, "/", L, "\r", appendLF = FALSE)
    # Observation l 
    data_l <- full_data[l, ]
    lat0 <- data_l$lat
    lon0 <- data_l$lon
    theta[l] <- data_l$x[[1]]
    # Remove observation l from full data
    reduced_data <- full_data[-l, ]
    # Fit the model without the l-th parameter 
    ok <- ordinaryKriging$new(reduced_data$x, reduced_data$lat, reduced_data$lon)
    ok$fit(model = cor_model)
    # Predict mean and variance of the l-th parameter
    theta_hat <- purrr::quietly(ok$predict)(lat = lat0, lon = lon0, n = n.neighb)$result
    m_theta[l] <- theta_hat$theta
    v_theta[l] <- theta_hat$var_theta
  }
  
  # CV score -> best minimize CVNLL 
  CVNLL <- sum(0.5 * log(2 * base::pi * v_theta) + (theta - m_theta)^2 / (2 * v_theta)) / L
  # CV-RMSE (simpler)
  CVRMSE <- sqrt(sum((theta - m_theta)^2) / L)
  
  # Results 
  dplyr::tibble(
    model = cor_model,
    formula = NA,
    CVNLL = CVNLL, 
    CVRMSE = CVRMSE, 
    theta = list(theta), 
    m_theta = list(m_theta),
    v_theta = list(v_theta)
  )
}


#' Leave-one-out cross-validation for k-nearest-neighbour kriging
#'
#' Performs leave-one-out cross-validation (LOOCV) for a k-nearest-neighbour
#' kriging interpolator. At each iteration, one observation is removed, the
#' interpolator is fitted on the remaining locations, and the omitted value
#' is predicted using the weighted average of its nearest neighbours.
#'
#' The function evaluates the interpolation accuracy using the
#' root mean squared prediction error (CVRMSE).
#'
#' @param x Numeric vector of observed values to interpolate.
#' @param lat Numeric vector of latitudes associated with \code{x}.
#' @param lon Numeric vector of longitudes associated with \code{x}.
#' @param n.neighb Integer specifying the number of nearest neighbours used
#'   in the interpolation. Default is \code{10}.
#' @param weight_fun Function used to compute neighbour weights as a function
#'   of distance. The function should take distances as input and return
#'   corresponding weights.
#' @param quiet Logical. If \code{FALSE}, prints the current LOOCV iteration.
#'   Default is \code{FALSE}.
#'
#' @details
#' For each location \eqn{\ell = 1, \dots, L}, the function:
#' \enumerate{
#'   \item removes observation \eqn{\ell};
#'   \item fits a \code{knnKriging} model on the remaining \eqn{L-1} locations;
#'   \item predicts the omitted value at
#'   \eqn{(\mathrm{lat}_\ell, \mathrm{lon}_\ell)};
#'   \item stores the predicted value.
#' }
#'
#' The cross-validation score is computed as
#' \deqn{
#' \mathrm{CVRMSE} =
#' \sqrt{\frac{1}{L} \sum_{\ell=1}^L
#' (\theta_\ell - \hat{\theta}_\ell)^2},
#' }
#' where \eqn{\theta_\ell} is the observed value and
#' \eqn{\hat{\theta}_\ell} is the leave-one-out prediction.
#'
#' @return A tibble with one row containing:
#' \describe{
#'   \item{\code{n.neighb}}{Number of neighbours used in interpolation.}
#'   \item{\code{weight_fun}}{List-column containing the weighting function.}
#'   \item{\code{CVRMSE}}{Cross-validation root mean squared error.}
#'   \item{\code{theta}}{List-column containing the observed values.}
#'   \item{\code{m_theta}}{List-column containing the LOOCV predicted values.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' lat <- runif(20, 44, 45)
#' lon <- runif(20, 11, 12)
#' x   <- rnorm(20)
#'
#' w_fun <- function(d) 1 / (d + 1e-6)
#'
#' cross_validate_knn_kriging(
#'   x = x,
#'   lat = lat,
#'   lon = lon,
#'   n.neighb = 8,
#'   weight_fun = w_fun
#' )
#' }
#'
#' @export
cross_validate_knn_kriging <- function(x, lat, lon, n.neighb = 10, weight_fun, quiet = FALSE){
  # Full dataset with all the locations
  full_data <- data.frame(x = x, lat = lat, lon = lon)
  L <- nrow(full_data)
  
  # True parameter 
  theta   <- vector("numeric", length = L)
  # Interpolated parameter 
  m_theta <- vector("numeric", length = L)
  for(l in 1:L) {
    if(!quiet) message(l, "/", L, "\r", appendLF = FALSE)
    # Observation l 
    data_l <- full_data[l, ]
    lat0 <- data_l$lat
    lon0 <- data_l$lon
    theta[l] <- data_l$x[[1]]
    # Remove observation l from full data
    reduced_data <- full_data[-l, ]
    # Fit the model without the l-th parameter 
    knn <- knnKriging$new(reduced_data$x, reduced_data$lat, reduced_data$lon)
    knn$fit(weight_fun = weight_fun, n.neighb = n.neighb)
    # Predict mean and variance of the l-th parameter
    m_theta[l] <- knn$predict(lat0, lon0, n = n.neighb)$theta
  }
  # CV-RMSE 
  CVRMSE <- sqrt(sum((theta - m_theta)^2) / L)
  # Results 
  dplyr::tibble(
    n.neighb = n.neighb,
    weight_fun = list(weight_fun),
    CVRMSE = CVRMSE, 
    theta = list(theta), 
    m_theta = list(m_theta)
  )
}

cross_validate_knn <- function(x, lat, lon, weight_fun, init.neighb = 10, init_theta = c(log(20000)), theta_params_fun = function(x) exp(x), max.neighb = 250){
  # Map constraint to uncontraint (n.neighb)
  theta_to_n.neighb <- function(max.neighb){
    force(max.neighb)
    function(theta){
      trunc(max.neighb * plogis(theta))
    }
  }
  theta_to_n.neighb <- theta_to_n.neighb(max.neighb)
  # Map unconstraint to contraint (n.neighb)
  n.neighb_to_theta <- function(max.neighb){
    force(max.neighb)
    function(n.neighb){
      qlogis(n.neighb/max.neighb)
    }
  }
  n.neighb_to_theta <- n.neighb_to_theta(max.neighb)
  # Loss function 
  loss_function <- function(params, weight_fun, theta_params_fun){
    # Number of neighborhoods
    n.neighb <- theta_to_n.neighb(params[1])
    # Range for weights
    theta_params <- theta_params_fun(params[-1])
    # Custom weight function 
    w_fun <- weight_fun(theta_params)
    # Perform a cross validation 
    cv <- cross_validate_knn_kriging(x, lat, lon,  n.neighb = n.neighb, w_fun, quiet = FALSE)
    # Loss CVRMSE
    print(paste0("n.neighb: ", n.neighb, " - Theta: ", theta_params,   " - CVRMSE: ", cv$CVRMSE))
    cv$CVRMSE
  }
  # ******************************************************************
  # Initial parameters
  init_params <- c(n.neighb_to_theta(init.neighb), init_theta)
  # Numerical optimization 
  opt <- optim(init_params, loss_function, weight_fun = weight_fun, theta_params_fun = theta_params_fun)
  # Optimal parameters
  params <- opt$par
  # Number of neighborhoods
  n.neighb <- theta_to_n.neighb(params[1])
  # Range for weights
  theta_params <- theta_params_fun(params[-1])
  # Custom weight function 
  w_fun <- weight_fun(theta_params)
  cat("------------------------------------------------------ \n")
  print(paste0("Best model params: n.neighb (", n.neighb, "), Theta (", paste0(theta_params, collapse = " - "),   ") - CVRMSE (", opt$value, ") \n"))
  cat("------------------------------------------------------ \n")
  # Evaluate the best model
  cv_best <- cross_validate_knn_kriging(x, lat, lon,  n.neighb = n.neighb, w_fun, quiet = FALSE)
  list(
    n.neighb = n.neighb, 
    theta_params = theta_params,
    init.neighb = init.neighb,
    init_theta = init_theta,
    theta_to_n.neighb = theta_to_n.neighb, 
    n.neighb_to_theta = n.neighb_to_theta, 
    theta_params_fun = theta_params_fun,
    cv = cv_best
  )
}





