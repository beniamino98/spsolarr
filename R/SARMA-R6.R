#' `SARMA` object
#'
#' @rdname SARMA
#' @name SARMA
#' @export
SARMA <- R6::R6Class("SARMA", 
                     public = list(
                       initialize = function(sar = TRUE, sma = FALSE){
                         # Store model's specification
                         private[["..sar"]] <- sar
                         private[["..sma"]] <- sma
                       },
                       fit = function(Y, W, dates, hsk = FALSE, monthly = FALSE){
                         private[["hsk"]] <- hsk
                         private[["monthly"]] <- monthly
                         # Fit the best parameters
                         sarma <- SARMA_fit(sar = private$..sar, sma = private$..sma, Y, W, dates, hsk = hsk, monthly = monthly)
                         # Store fitted parameters 
                         private[["..rho"]] <- sarma$sar$rho
                         private[["..lambda"]] <- sarma$sma$lambda
                         private[["..Sigma2"]] <- sarma$Sigma2
                         private[["..log.likelihoods"]] <- sarma$log.likelihoods
                       },
                       filter_to_Z = function(Y, W, dates){
                         # Number of observations per locations
                         N <- nrow(Y)
                         # Number of locations 
                         L <- ncol(Y)
                         # Product
                         iBA <- self$iBA(W)
                         # Fitted esiduals
                         eps_hat <-  Y %*% t(iBA)
                         
                         # Standardized residuals
                         z_hat <- eps_hat 
                         if (private$monthly){
                           nmonths <- lubridate::month(dates)
                           for(nmonth in unique(nmonths)){
                             Sigma2_m <- self$Sigma2[[nmonth]]
                             idx_month <- nmonths == nmonth
                             # Standardized residuals 
                             z_hat[idx_month,] <- eps_hat[idx_month,] %*% diag(1/sqrt(diag(Sigma2_m)))
                           }
                         } else {
                           z_hat <- eps_hat %*% diag(1/sqrt(diag(self$Sigma2)))
                         }
                         z_hat
                       },
                       filter_to_Y = function(Z, W, dates){
                         # Number of observations per locations
                         N <- nrow(Z)
                         # Number of locations 
                         L <- ncol(Z)
                         # From z_hat to eps_hat
                         eps_hat <- Z
                         if (private$monthly) {
                           nmonths <- lubridate::month(dates)
                           for(nmonth in unique(nmonths)){
                             Sigma2_m <- self$Sigma2[[nmonth]]
                             idx_month <- nmonths == nmonth
                             # Standardized residuals 
                             eps_hat[idx_month,] <- Z[idx_month,] %*% diag(sqrt(diag(Sigma2_m)))
                           }
                         } else {
                           eps_hat <- Z %*% diag(sqrt(diag(self$Sigma2)))
                         }
                         # Product
                         BiA <- self$BiA(W)
                         # Fitted values Y
                         eps_hat %*% t(BiA)
                       },
                       iBA = function(W){
                         L <- nrow(W)
                         # Identity matrix 
                         I <- diag(1, L, L)
                         # SAR matrix 
                         A <- I - self$rho * W 
                         # SMA matrix 
                         B <- I - self$lambda * W 
                         # Product
                         solve(B) %*% A
                       },
                       BiA = function(W){
                         L <- nrow(W)
                         # Identity matrix 
                         I <- diag(1, L, L)
                         # SAR matrix 
                         A <- I - self$rho * W 
                         # SMA matrix 
                         B <- I - self$lambda * W 
                         # Product
                         solve(A) %*% B
                       }
                     ), 
                     private = list(
                       ..sar = TRUE,
                       ..sma = TRUE, 
                       ..rho = 0,
                       ..lambda = 0,
                       ..Sigma2 = list(),
                       ..log.likelihoods = NA,
                       
                       hsk = FALSE, 
                       monthly = FALSE
                     ), 
                     active = list(
                       rho = function(){
                         private$..rho
                       },
                       lambda = function(){
                         private$..lambda
                       },
                       loglik = function(){
                         sum(private$..log.likelihoods)
                       },
                       Sigma2 = function(){
                         private$..Sigma2
                       }
                     ))
