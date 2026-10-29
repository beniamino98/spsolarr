#' `simpleKriging` object
#'
#' @rdname simpleKriging
#' @name simpleKriging
#' @export
simpleKriging <- R6::R6Class(
  "simpleKriging",
  public = list(
    #' @field var.name Variable name 
    var.name = "theta",
    #' @field model.name Reference model name 
    model.name = "Sph",
    #' @field formula Reference formula for the linear trend 
    formula = "x ~ 1",
    #' @description
    #' Initialization 
    #' @param x numeric vector to fit
    #' @param lat numeric vector of latitudes in degrees
    #' @param lon numeric vector of longitudes in degrees
    #' @param var.name character, custom parameter label 
    initialize = function(x, lat, lon, var.name = "theta") {
      # Project coordinates 
      spdf <- project_CRS(lat, lon)
      # Add variable name 
      spdf@data[[var.name]]  <- x
      spdf@data[["x"]] <- x
      # Store data 
      private$..data <- spdf
      # Store variable name 
      self$var.name <- var.name
    },
    #' @description
    #' Fit a kriging model 
    #' @param mu_theta numeric scalar, fixed mean of the parameter. 
    #' @param model character, variogram model. 
    #' @param gls gstat args
    #' @param ... extra gstat args. 
    fit = function(mu_theta = 0, model = "Sph", gls = FALSE, ...){
      # Variable name 
      data <- self$data
      # Reference model 
      self$model.name <- model
      # Add mean of mu_theta 
      data@data[["mu_theta"]] <- mu_theta
      # Remove mean from the variable 
      data@data[["x"]] <- data@data$x - mu_theta
      # Empirical variogram 
      vg_emp <- gstat::variogram(stats::as.formula(self$formula), data = data)
      # Fitted variogram 
      vg_model <- gstat::fit.variogram(vg_emp, model = gstat::vgm(model = model, nugget = 1e-6), fit.kappa = TRUE)
      # G-stat object for prediction
      # Simple kriging uses beta = known mean (after drift removal, if any)
      g <- gstat::gstat(
        id    = "x",
        formula = stats::as.formula(self$formula),
        data  = data,
        model = vg_model,
        beta  = mu_theta,
        set   = if (isTRUE(gls)) list(gls = 1) else list()
      )
      private$..mu_theta <- mu_theta
      private$..vg_emp   <- vg_emp
      private$..vg_model <- vg_model
      private$..g        <- g
    },
    #' @description
    #' Predict 
    #' @param lat numeric vector of latitudes in degrees
    #' @param lon numeric vector of longitudes in degrees
    #' @param n numeric, number of neighbouroods.
    predict = function(lat, lon, n = 4){
      # Safety check 
      stopifnot(!is.null(private$..g))
      # Projected coordnates 
      new_sp <- project_CRS(lat, lon)
      # gstat wants a data slot
      new_sp@data[["x"]] <- NA_real_
      # Local kriging via nmax
      pr <- gstat:::predict.gstat(private$..g, newdata = new_sp, nmax = n)
      # Predicted values 
      out <- dplyr::tibble(lat = lat, lon = lon, interpolated = TRUE)
      # Assign a name equal to the variable name 
      out[[self$var.name]] <- pr@data[["x.pred"]] + private$..mu_theta
      out[[paste0("var_", self$var.name)]] <- pr@data[["x.var"]]
      out
    },
    #' @description
    #' Print method for the R6 class `simpleKriging`. 
    print = function(){
      cat("------------------- simpleKriging ------------------- \n")
      cat(paste0("Variogram model: ", self$model.name), "\n")
      cat(paste0("Variogram formula: ", self$formula), "\n")
      cat(paste0("Mu theta: ", private$..mu_theta), "\n")
    }
  ),
  private = list(
    ..data = NULL,
    ..vg_emp = NULL,
    ..vg_model = NULL,
    ..g = NULL,
    ..mu_theta = 0
  ),
  active = list(
    #' @field data parameters data. 
    data = function() private$..data,
    #' @field vg_model variogram model. 
    vg_model = function() private$..vg_model,
    #' @field vg_emp empiric variogram. 
    vg_emp = function() private$..vg_emp,
    #' @field g gstat object
    g = function() private$..g
  )
)


#' `ordinaryKriging` object
#'
#' @rdname ordinaryKriging
#' @name ordinaryKriging
#' @export
ordinaryKriging <- R6::R6Class(
  "ordinaryKriging",
  inherit = simpleKriging,
  public = list(
    #' @description
    #' Fit a kriging model 
    #' @param model character, variogram model. 
    #' @param gls gstat args
    #' @param ... extra gstat args. 
    fit = function(model = "Sph", gls = FALSE, ...) {
      # Reference model 
      self$model.name <- model
      # In OK: unknown constant mean, handled only by krige formula x ~ 1
      super$fit(mu_theta = 0, model = model, gls = gls)
      vg_emp <- private$..vg_emp
      vg_model <- private$..vg_model
      # Ensure non-zero psill 
      if (vg_model$psill[1] == 0){
        vg_model$psill[1] <- 10e-6
      }
      # BUT: in OK must NOT pass beta. So rebuild g without beta:
      data <- self$data
      private$..g <- gstat::gstat(
        id = "x",
        formula = as.formula(self$formula),
        data = data,
        model = vg_model,
        set = if (isTRUE(gls)) list(gls = 1) else list()
      )
    },
    #' @description
    #' Print method for the R6 class `ordinaryKriging`. 
    print = function(){
      cat("------------------- ordinaryKriging ------------------- \n")
      cat(paste0("Variogram model: ", self$model.name), "\n")
      cat(paste0("Variogram formula: ", self$formula), "\n")
    }
  )
)

#' `universalKriging` object
#'
#' @rdname universalKriging
#' @name universalKriging
#' @export
universalKriging <- R6::R6Class(
  "universalKriging",
  inherit = simpleKriging,
  public = list(
    #' @description
    #' Fit the kriging model 
    #' @param formula character formula 
    #' @param model character, variogram model. 
    #' @param gls gstat args
    #' @param ... extra gstat args. 
    fit = function(formula = "x ~ 1 + lat + lon", model = "Sph", gls = FALSE, ...) {
      # Data 
      data <- self$data
      # Reference model 
      self$model.name <- model
      # Reference model 
      self$formula <- stats::as.formula(formula)
      # add coordinate covariates for drift (projected)
      xy <- sp::coordinates(data)
      data@data[["lon"]] <- xy[,1]
      data@data[["lat"]] <- xy[,2]
      # Empiric variogram 
      vg_emp <- gstat::variogram(self$formula, data = data)
      # Fitted variogram 
      vg_model <- gstat::fit.variogram(vg_emp, model = gstat::vgm(model = model, nugget = 1e-6), fit.kappa = TRUE)
      if (vg_model$psill[1] == 0){
        vg_model$psill[1] <- 10e-6
      }
      # gstat object 
      g <- gstat::gstat(
        id = "x",
        formula = self$formula,
        data = data,
        model = vg_model,
        set = if (isTRUE(gls)) list(gls = 1) else list()
      )
      private$..vg_emp   <- vg_emp
      private$..vg_model <- vg_model
      private$..g        <- g
    },
    #' @description
    #' Predict 
    #' @param lat numeric vector of latitudes in degrees
    #' @param lon numeric vector of longitudes in degrees
    #' @param n numeric, number of neighborhoods.
    predict = function(lat, lon, n = 4) {
      stopifnot(!is.null(private$..g))
      # Project coordinates 
      new_sp <- project_CRS(lat, lon)
      # Add covariates
      xy <- sp::coordinates(new_sp)
      new_sp@data[["lon"]] <- xy[,1]
      new_sp@data[["lat"]] <- xy[,2]
      # Add target data 
      new_sp@data[["x"]] <- NA_real_
      # Predict 
      pr <- gstat:::predict.gstat(private$..g, newdata = new_sp, nmax = n, quiet = TRUE)
      out <- dplyr::tibble(lat = lat, lon = lon, interpolated = TRUE)
      out[[self$var.name]] <- pr@data[["x.pred"]]
      out[[paste0("var_", self$var.name)]] <- pr@data[["x.var"]]
      out
    },
    #' @description
    #' Print method for the R6 class `universalKriging`. 
    print = function(){
      cat("------------------- universalKriging ------------------- \n")
      cat(paste0("Variogram model: ", self$model.name), "\n")
      cat(paste0("Variogram formula: ", self$formula), "\n")
    }
  )
)

#' `knnKriging` object
#'
#' kNN / IDW-style interpolation of a spatial parameter using neighborhood values only.
#' This is a baseline competitor to kriging and is useful also as a robust fallback.
#'
#' @rdname knnKriging
#' @name knnKriging
#' @export
knnKriging <- R6::R6Class("knnKriging",
                          public = list(
                            #' @field var.name Variable name 
                            var.name = "theta",
                            #' @field model.name Reference model name 
                            model.name = "knn",
                            #' @field formula Reference formula for the linear trend 
                            formula = NA_character_,
                            #' @description
                            #' Initialization 
                            #' @param x numeric vector, parameter values at observed locations
                            #' @param lat numeric vector of latitudes in degrees
                            #' @param lon numeric vector of longitudes in degrees
                            #' @param var.name character, custom parameter label 
                            initialize = function(x, lat, lon, var.name = "theta"){
                              data_sp <- project_CRS(lat, lon)
                              # Store variable name 
                              self$var.name <- var.name
                              data_sp[[var.name]] <- x
                              private[["..data"]] <- data_sp
                            },
                            #' @description 
                            #' default neighborhoods method (same logic as your kriging classes)
                            #' @param lat numeric vector of latitudes in degrees
                            #' @param lon numeric vector of longitudes in degrees
                            #' @param n numeric, number of neighborhoods.
                            neighborhoods = function(lat, lon, n = 4){
                              coords <- dplyr::tibble(lat = lat, lon = lon)
                              coords_proj <- project_CRS(coords$lat, coords$lon)
                              
                              data <- self$data
                              res <- vector("list", length(lat))
                              isKnown <- logical(length(lat))
                              
                              for (i in seq_along(lat)) {
                                coords_k <- coords_proj[i, ]
                                h <- sp::spDistsN1(data, coords_k, longlat = FALSE)
                                idx_ngb <- order(h, decreasing = FALSE)[1:n]
                                
                                res[[i]] <- data[idx_ngb, ]
                                res[[i]]$h <- h[idx_ngb]
                                isKnown[i] <- any(res[[i]]$h == 0)
                              }
                              
                              coords_out <- dplyr::as_tibble(sp::coordinates(coords_proj))
                              coords_out$isKnown <- isKnown
                              coords_out$ngb <- res
                              coords_out
                            },
                            #' @description
                            #' Store the weight function / its parameters. No variogram fit here.
                            #' @param weight_fun function(h, ngb, lat0, lon0, ...) -> numeric weights (length n)
                            #' @param ... optional extra args forwarded to weight_fun
                            fit = function(weight_fun = NULL, ...){
                              if (is.null(weight_fun)) {
                                # Default: IDW with power p=2 and small eps for stability
                                weight_fun <- function(beta = 2) {
                                  force(beta); 
                                  function(h){
                                    1 / (h)^beta
                                  }
                                }
                                private[["..weight_fun"]] <- weight_fun(2)
                              } else {
                                private[["..weight_fun"]] <- weight_fun
                              }
                            },
                            #' @description
                            #' Predict by weighted average of neighbor parameter values.
                            #' @param lat numeric vector of latitudes in degrees
                            #' @param lon numeric vector of longitudes in degrees
                            #' @param n numeric, number of neighbouroods.
                            #' @param coords_ngb optional. 
                            predict = function(lat, lon, n = 4, coords_ngb){
                              var.name <- self$var.name
                              
                              if (missing(coords_ngb)) {
                                coords_ngb <- self$neighborhoods(lat, lon, n)
                              }
                              
                              preds <- numeric(nrow(coords_ngb))
                              
                              for (i in seq_len(nrow(coords_ngb))) {
                                coords_ngb_i <- coords_ngb[i, ]
                                # exact copy if location coincides with known point
                                if (isTRUE(coords_ngb_i$isKnown)) {
                                  preds[i] <- coords_ngb_i$ngb[[1]][[var.name]][1]
                                  next
                                }
                                
                                ngb <- coords_ngb_i$ngb[[1]]
                                theta_ngb <- ngb[[var.name]]
                                
                                # user-defined weights
                                w <- private$..weight_fun(ngb$h)
                                w <- w / sum(w)
                                # safety checks
                                w[!is.finite(w)] <- 0
                                if (length(w) != length(theta_ngb)) {
                                  stop("weight_fun must return a vector of length equal to number of neighbors.")
                                }
                                s <- sum(w)
                                if (s <= 0) {
                                  # fallback: uniform weights
                                  w <- rep(1/length(theta_ngb), length(theta_ngb))
                                } else {
                                  w <- w / s
                                }
                                preds[i] <- sum(w * theta_ngb)
                              }
                              
                              out <- dplyr::tibble(lat = lat, lon = lon, interpolated = !coords_ngb$isKnown)
                              out[[var.name]] <- preds
                              out[[paste0("var_", var.name)]] <- 0
                              out
                            },
                            #' @description
                            #' Print method for the R6 class `knnKriging`. 
                            print = function(){
                              cat("------------------- knnInterpolation ------------------- \n")
                              cat(paste0("Variogram model: ", self$model.name, "\n"))
                              cat(paste0("Variogram formula: ", self$formula, "\n"))
                              cat(paste0("Weighting function: ", as.character(body(private$..weight_fun)), "\n"))
                            }
                          ),
                          private = list(
                            ..data = NA,
                            ..weight_fun = NULL,
                            ..weight_args = list()
                          ),
                          active = list(
                            #' @field data parameters data. 
                            data = function() private$..data
                            )
)


