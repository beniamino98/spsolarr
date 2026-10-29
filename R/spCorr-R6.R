#' `spCorr` object
#'
#' @rdname spCorr
#' @name spCorr
#' @export
spCorr <- R6::R6Class("spCorr", 
                      public = list(
                        #' @description
                        #' Initialization 
                        #' @param model.name character, model name
                        initialize = function(model.name = "exp"){
                          private$..params <- sp_cor_params(model.name)$params
                          private$..model.name <- model.name
                        },
                        #' @description
                        #' Compute correlation from a distance
                        #' @param h distance
                        rho = function(h){
                          # Reference model
                          model.name <- private$..model.name
                          # Parameters
                          params = private$..params
                          # Function name 
                          fun.name <- paste0("sp_cor_", model.name)
                          # Rho 
                          do.call(fun.name, args = list(h = h, params = params))
                        },
                        #' @description
                        #' Update the parameters
                        #' @param coefficients  new parameters 
                        update = function(coefficients){
                          private$..params <- coefficients
                        },
                        #' @description
                        #' Print method for the R6 class `spCorr`. 
                        print = function(){
                          params <- paste0(names(private$..params), ": ", private$..params, collapse = "\n ")
                          cat("-------------- spatial Correlation -------------- \n")
                          cat("Model: ", private$..model.name, "\n")
                          
                          cat("Parameters: \n", params, "\n")
                          cat("------------------------------------------------- \n")
                        }
                      ), 
                      private = list(
                        ..model.name = "exp",
                        ..params = c(range = 1)
                      ), 
                      active = list(
                        #' @field params params for correlation 
                        params = function(){
                          private$..params
                        },
                        #' @field model.name  model name
                        model.name = function(){
                          private$..model.name
                        },
                        #' @field params_star tranasformed parameters 
                        params_star = function(){
                          sp_cor_params_to_zeta(private$..params, model.name = private$..model.name)
                        }
                      )
)

