#' `SARMA_contingency_mat` 
#'
#' @rdname SARMA_contingency_mat
#' @name SARMA_contingency_mat
#' @export
SARMA_contingency_mat <- function(D, n.neighb = 4, dist_fun = function(x)1/x, 
                                  spectral.normalize = TRUE){
  # Number of points
  L <- nrow(D)
  # Initialization 
  W <- D
  for(l in 1:L){
    # Detect closest neightbouroods 
    idx_neigh <- order(D[l,], decreasing = FALSE)[2:(n.neighb+1)]
    W[l,idx_neigh] <- dist_fun(W[l,idx_neigh])
    W[l, -idx_neigh] <- 0
  }
  # Spectral normalization 
  if (spectral.normalize){
    W_lam <- eigen(W)$values
    lam_max <- max(abs(W_lam))
    W <- W / lam_max
  } else {
    # Row normalization 
    W <- W / rowSums(W)
  }
  return(W)
}

#' `SARMA_filter` 
#'
#' @rdname SARMA_filter
#' @name SARMA_filter
#' @export
SARMA_filter <- function(rho, lambda, Y, W, dates, hsk = FALSE, monthly = FALSE){
  # Number of observations per locations
  N <- nrow(Y)
  # Number of locations 
  L <- ncol(Y)
  # Identity matrix 
  I <- diag(1, L, L)
  # SAR matrix 
  A <- I - rho * W 
  # SMA matrix 
  B <- I - lambda * W 
  # Product
  iBA <- solve(B) %*% A
  # Fitted esiduals
  eps_hat <-  Y %*% t(iBA)

  # Initialize standardized residuals
  z_hat <- eps_hat 
  # Estimation of the variance 
  if (monthly){
    nmonths <- lubridate::month(dates)
    Sigma2 <- list()
    for(nmonth in unique(nmonths)){
      Sigma2_m <- matrix(0, L, L)
      idx_month <- nmonths == nmonth
      # Heteroskedasticity 
      if (hsk) {
        diag(Sigma2_m) <- apply(eps_hat[idx_month,], 2, function(x) mean(x^2))
      } else {
        diag(Sigma2_m) <- sum(apply(eps_hat[idx_month,], 2, function(x) sum(x^2))) / (sum(idx_month)*L)
      }
      Sigma2[[nmonth]] <- Sigma2_m
      # Standardized residuals 
      z_hat[idx_month,] <- eps_hat[idx_month,] %*% diag(1/sqrt(diag(Sigma2_m)))
    }
    names(Sigma2) <- as.character(lubridate::month(unique(nmonths), label = TRUE))
  } else {
    Sigma2 <- matrix(0, L, L)
    if (hsk) {
      diag(Sigma2) <- apply(eps_hat, 2, function(x) mean(x^2))
    } else {
      diag(Sigma2) <- sum(apply(eps_hat, 2, function(x) sum(x^2))) / (N*L)
    }
    # Standardized residuals 
    z_hat <- eps_hat %*% diag(1/sqrt(diag(Sigma2)))
  }
  
  # Log-likelihood 
  iiBA <- solve(iBA)
  if (monthly){
    log.likelihoods <- numeric(N)
    for(nmonth in unique(nmonths)){
      Sigma2_m <- Sigma2[[nmonth]]
      idx_month <- nmonths == nmonth
      # Variance-covariance matrix 
      V <- t(iiBA) %*% Sigma2_m %*% iiBA
      log.likelihoods[idx_month] <- mvtnorm::dmvnorm(Y[idx_month,], sigma = V, log = TRUE)
    }
  } else {
    # Variance-covariance matrix 
    V <- t(iiBA) %*% Sigma2 %*% iiBA
    # Log-likelihoods
    log.likelihoods <- mvtnorm::dmvnorm(Y, sigma = V, log = TRUE)
  }
  
  structure(
    list(
      Y = Y, 
      W = W, 
      eps_hat = eps_hat, 
      z_hat = z_hat,
      sar = list(A = A, rho = rho),
      sma = list(B = B, lambda = lambda),
      Sigma2 = Sigma2,
      log.likelihoods = log.likelihoods,
      dims = list(N = N, L = L),
      hsk = hsk,
      monthly = monthly
    )
  )
}

#' `SARMA_fit` 
#'
#' @rdname SARMA_fit
#' @name SARMA_fit
#' @export
SARMA_fit <- function(sar = TRUE, sma = TRUE, Y, W, dates, hsk = TRUE, monthly = FALSE){
  # Loss function 
  loss_function <- function(params){
    if (sar == TRUE & sma == TRUE) {
      sarma <- SARMA_filter(params[1], params[2], Y, W, dates, hsk = hsk, monthly = monthly)
    } else {
      sarma <- SARMA_filter(params[1], 0, Y, W, dates, hsk = hsk, monthly = monthly)
    }
    # Log-likelihood 
    loglik <- sum(sarma$log.likelihoods)
    # Verbose message
    print(paste0("Log-lik: ", loglik, " Params: ", params))
    return(-loglik)
  }
  if (sar == TRUE & sma == TRUE) {
    opt <- optim(c(0, 0), loss_function, 
                 method = "L-BFGS-B", lower = c(0, -0.99), upper = c(0.99, 0.99))
  } else {
    opt <- optim(0, loss_function, 
                 method = "Brent", lower = c(0), upper = c(0.99))
  }
  # Filter 
  params <- opt$par
  if (sar == TRUE & sma == TRUE) {
    sarma <- SARMA_filter(params[1], params[2], Y, W, dates, hsk = hsk, monthly = monthly)
  } else {
    sarma <- SARMA_filter(params[1], 0, Y, W, dates, hsk = hsk, monthly = monthly)
  }
  return(sarma)
}

#' `SARMA_test_ACF` 
#'
#' @rdname SARMA_test_ACF
#' @name SARMA_test_ACF
#' @export
SARMA_test_ACF <- function(cr_z_hat, Tn){
  # Number of locations
  L <- ncol(cr_z_hat)
  # Number of observations
  Tn <- Tn
  # Leave only unique correlations
  cr_z <- cr_z_hat
  cr_z[upper.tri(cr_z)] <- 0
  diag(cr_z) <- 0
  # Statistic test
  stat <- sqrt(2/(L*(L-1))) * sum(cr_z) / sqrt(Tn)
  # P.value
  p.value <- 1-pnorm(abs(stat))
  dplyr::tibble(method = "Pesaran CD test", stat = stat, p.value = p.value)
}


