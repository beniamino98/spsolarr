#' `sp_cor` Models 
#'
#' @rdname sp_cor
#' @name sp_cor
#' @aliases sp_cor_av.models
#' @aliases sp_cor_params
#' @aliases sp_cor_params_to_zeta
#' @aliases sp_cor_params_to_phi
#' @export
sp_cor_av.models <- function(){
  c("exp", "gau", "powExp", "sph", "ratQuad", "matern", "cauchy", "dagum", "wendland", "mixExp", "mixCauchy")
}
#' @rdname sp_cor
#' @export
sp_cor_params <- function(model.name = "Exp"){
  # Base parameters
  params <- c(range = 10000)
  # Available models
  model.name <- match.arg(model.name, choices = sp_cor_av.models())
  # Extra parameters
  if (model.name == "ratQuad"){
    params <- c(params, alpha = 0.5)
  } else if (model.name == "matern"){
    params <- c(params, nu = 0.5)
  } else if (model.name == "powExp"){
    params <- c(params, kappa = 0.5)
  } else if (model.name == "cauchy"){
    params <- c(params, alpha = 0.5, beta = 0.5)
  } else if (model.name == "dagum"){
    params <- c(params, alpha = 0.5, beta = 0.5)
  } else if (model.name == "mixExp"){
    params <- c(params, range2 = 10000, weight = 0.5)
  } else if (model.name == "mixCauchy"){
    params <- c(params, range2 = 10000, weight = 0.5, alpha = 1.2, beta = 1.2, alpha2 = 1.2, beta2 = 1.2)
  } 
  list(
    model.name = model.name,
    params = params
  )
}

#' @rdname sp_cor
#' @export
sp_cor_params_to_zeta <- function(params, model.name = "Exp"){
  # Available models
  model.name <- match.arg(model.name, choices = sp_cor_av.models())
  params_star <- c(range_star = log(params[["range"]]))
  # Extra parameters
  if (model.name == "ratQuad"){
    # Alpha > 0
    params_star <- c(params_star, alpha_star = log(params[["alpha"]]))
  } else if (model.name == "matern"){
    # Nu > 0
    params_star <- c(params_star, nu_star = log(params[["nu"]]))
  } else if (model.name == "powExp"){
    # 0 < kappa <= 2
    params_star <- c(params_star, kappa_star = qlogis(params[["kappa"]]/2))
  } else if (model.name == "cauchy"){
    # 0 < alpha <= 2, beta > 0
    params_star <- c(params_star, alpha_star = qlogis(params[["alpha"]]/2), beta_star = log(params[["beta"]]))
  } else if (model.name == "dagum"){
    # alpha > 0, beta > 0
    params_star <- c(params_star, alpha_star = log(params[["alpha"]]), beta_star = log(params[["beta"]]))
  } else if (model.name == "mixExp"){
    params_star <- c(params_star, range2_star = log(params[["range2"]]), weight_star = qlogis(params[["weight"]]))
  } else if (model.name == "mixCauchy"){
    params_star <- c(params_star, range2_star = log(params[["range2"]]), weight_star = qlogis(params[["weight"]]), 
                     alpha_star = qlogis(params[["alpha"]]/2), beta_star = log(params[["beta"]]), 
                     alpha2_star = qlogis(params[["alpha2"]]/2), beta2_star = log(params[["beta2"]]))
  } 
  list(
    model.name = model.name,
    params_star = params_star
  )
}

#' @rdname sp_cor
#' @export
sp_cor_params_to_phi <- function(params_star, model.name = "Exp"){
  # Available models
  model.name <- match.arg(model.name, choices = sp_cor_av.models())
  params <- c(range = exp(params_star[["range_star"]]))
  # Extra parameters
  if (model.name == "ratQuad"){
    # Alpha > 0
    params <- c(params, alpha = exp(params_star[["alpha_star"]]))
  } else if (model.name == "matern"){
    # Nu > 0
    params <- c(params, nu = exp(params_star[["nu_star"]]))
  } else if (model.name == "powExp"){
    # 0 < kappa <= 2
    params <- c(params, kappa = 2*plogis(params_star[["kappa_star"]]))
  } else if (model.name == "cauchy"){
    # 0 < alpha <= 2, beta > 0
    params <- c(params, alpha = 2*plogis(params_star[["alpha_star"]]), beta = exp(params_star[["beta_star"]]))
  } else if (model.name == "dagum"){
    # alpha > 0, beta > 0
    params <- c(params, alpha = exp(params_star[["alpha_star"]]), beta = exp(params_star[["beta_star"]]))
  } else if (model.name == "mixExp"){
    params <- c(params, range2 = exp(params_star[["range2_star"]]), weight = plogis(params_star[["weight_star"]]))
  } else if (model.name == "mixCauchy"){
    params <- c(params, range2 = exp(params_star[["range2_star"]]), weight = plogis(params_star[["weight_star"]]),
                alpha = 2*plogis(params_star[["alpha_star"]]), beta = exp(params_star[["beta_star"]]), 
                alpha2 = 2*plogis(params_star[["alpha2_star"]]), beta2 = exp(params_star[["beta2_star"]]))
  } 
  list(
    model.name = model.name,
    params = params
  )
}

#' Exponential spatial correlation
#'
#' @param h Matrix of vector of distances.
#' @param params Numeric named vector of parameters.
#' @note Version 1.0.0.
#' @details Implement the functions:
#' \deqn{\rho_{\exp}(h, \phi^{\exp}) = \exp\left(\frac{h}{\phi^{\exp}}\right)}
#' 
#' @rdname sp_cor_exp
#' @name sp_cor_exp
#' @export
sp_cor_exp <- function(h, params){
  range <- params["range"]
  h_range <- h / range
  exp(-h_range)
}

#' Gaussian spatial correlation
#'
#' @param h Matrix of vector of distances.
#' @param params Numeric named vector of parameters.
#' @note Version 1.0.0.
#' @details Implement the functions:
#' \deqn{\rho_{\text{gau}}(h, \phi^{\exp}) = \exp\left\{\left(\frac{h}{\phi^{\text{gau}}}\right)^2 \right\}}
#' 
#' @rdname sp_cor_gau
#' @name sp_cor_gau
#' @export
sp_cor_gau <- function(h, params){
  range <- params["range"]
  h_range <- h / range
  exp(-h_range^2)
}

#' Powered exponential spatial correlation
#'
#' @param h Matrix of vector of distances.
#' @param params Numeric named vector of parameters.
#' @note Version 1.0.0.
#' 
#' @rdname sp_cor_powExp
#' @name sp_cor_powExp
#' @export
sp_cor_powExp <- function(h, params){
  range <- params["range"]
  kappa <- params["kappa"]
  h_range <- h / range
  exp(-h_range^kappa)
}

#' Spherical spatial correlation
#'
#' @param h Matrix of vector of distances.
#' @param params Numeric named vector of parameters.
#' @note Version 1.0.0.
#' @details Implement the functions:
#' \deqn{\rho_{\text{sph}}(h, \phi^{\text{sph}}) =  1 - \frac{3}{2}\frac{h}{\phi^{\text{sph}}} + \frac{1}{2}\left(\frac{h}{\phi^{\text{sph}}}\right)^3}
#' 
#' @rdname sp_cor_sph
#' @name sp_cor_sph
#' @export
sp_cor_sph <- function(h, params){
  range <- params["range"]
  h_range <- h / range
  rho <- 1 - 1.5 * h_range + 0.5*h_range^3
  rho[h > range] <- 0
  rho
}

#' Rational quadratic spatial correlation
#'
#' @param h Matrix of vector of distances.
#' @param params Numeric named vector of parameters.
#' @note Version 1.0.0.
#' 
#' @rdname sp_cor_ratQuad
#' @name sp_cor_ratQuad
#' @export
sp_cor_ratQuad <- function(h, params){
  range <- params["range"]
  alpha <- params["alpha"]
  h_range <- h^2 / (alpha * range^2)
  (1 + h_range)^(-alpha)
}

#' Matern spatial correlation
#'
#' @param h Matrix of vector of distances.
#' @param params Numeric named vector of parameters.
#' @note Version 1.0.0.
#' 
#' @rdname sp_cor_matern
#' @name sp_cor_matern
#' @export
sp_cor_matern <- function(h, params){
  range <- params["range"]
  nu <- params["nu"]
  h_range <- h / range
  rho <- 1 / (2^(nu-1) * base::gamma(nu)) * (h_range^nu) * besselK(h_range, nu)
  rho[h == 0] <- 1
  rho
}

#' Generalized Cauchy spatial correlation
#'
#' @param h Matrix of vector of distances.
#' @param params Numeric named vector of parameters.
#' @note Version 1.0.0.
#' 
#' @rdname sp_cor_cauchy
#' @name sp_cor_cauchy
#' @export
sp_cor_cauchy <- function(h, params){
  range <- params["range"]
  alpha <- params["alpha"]
  beta <- params["beta"]
  h_range <- h / range
  (1 + h_range^alpha)^(-beta/alpha)
}

#' Generalized Cauchy spatial correlation
#'
#' @param h Matrix of vector of distances.
#' @param params Numeric named vector of parameters.
#' @note Version 1.0.0.
#' 
#' @rdname sp_cor_dagum
#' @name sp_cor_dagum
#' @export
sp_cor_dagum <- function(h, params){
  range <- params["range"]
  alpha <- params["alpha"]
  beta <- params["beta"]
  h_range <- h / range
  (1 - (1 + h_range^alpha)^(-beta))
}

#' Generalized Wenland spatial correlation
#'
#' @param h Matrix of vector of distances.
#' @param params Numeric named vector of parameters.
#' @note Version 1.0.0.
#' 
#' @rdname sp_cor_wendland 
#' @name sp_cor_wendland
#' @export
sp_cor_wendland <- function(h, params){
  range <- params["range"]
  h_range <- h / range
  rho <- (1 - h_range)^4 * (4 * h_range + 1)
  rho[h > range | 1 - h_range < 0] <- 0
  rho
}

#' Exponential spatial correlation
#'
#' @param h Matrix of vector of distances.
#' @param params Numeric named vector of parameters.
#' @note Version 1.0.0.
#' @details Implement the functions:
#' \deqn{\rho_{\exp}(h, \phi^{\exp}) = \exp\left(\frac{h}{\phi^{\exp}}\right)}
#' 
#' @rdname sp_cor_mixExp
#' @name sp_cor_mixExp
#' @export
sp_cor_mixExp <- function(h, params){
  range <- params["range"]
  range2 <- params["range2"]
  weight <- params["weight"]
  h_range <- h / range
  h_range2 <- h / range2
  weight * exp(-h_range) + (1 - weight) * exp(-h_range2)
}

#' Mixed cauchy 
#'
#' @param h Matrix of vector of distances.
#' @param params Numeric named vector of parameters.
#' @note Version 1.0.0.
#' @rdname sp_cor_mixCauchy
#' @name sp_cor_mixCauchy
#' @export
sp_cor_mixCauchy <- function(h, params){
  range <- params["range"]
  alpha <- params["alpha"]
  beta <- params["beta"]
  weight <- params["weight"]
  range2 <- params["range2"]
  alpha2 <- params["alpha2"]
  beta2 <- params["beta2"]
  h_range <- h / range
  h_range2 <- h / range2
  
  weight * (1 + h_range^alpha)^(-beta/alpha) + (1 - weight) * (1 + h_range2^alpha2)^(-beta2/alpha2)
}








