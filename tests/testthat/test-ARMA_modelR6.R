# ***********************************************************************************************
#                                  Test for class "ARMA_modelR6"
# ***********************************************************************************************
# Time series to use for fitting
x <- rnorm(1000)

# Initialize the model
arma <- ARMA_modelR6$new(maOrder = 1, arOrder = 1)
# Fit the model
arma$fit(x)

# Fitted residuals from arima
eps_arima <- arma$model$residuals
# Fitted residuals with ARMA$filter
eps_hat <- x - arma$filter(x)
sum(abs(eps_hat - eps_arima))

#private <- arma$.__enclos_env__$private
#self <- arma$.__enclos_env__$self

# arima model
arma$model

# Slots
## Orders
arma$arOrder
arma$maOrder
arma$order
## Parameters
arma$intercept
arma$phi
arma$theta
arma$coefficients
## Parameters std. errors
arma$std.errors
# Tidy ouptut
arma$tidy
## Matrices
arma$A
arma$b
## Estimated variance of the residuals
arma$sigma2

# Functions
## Expectation
arma$expectation(1, X0 = c(0.2, 0.1))
arma$expectation(10, X0 = c(0.2, 0.1))
arma$expectation(30, X0 = c(0.2, 0.1))
## Variance
arma$variance(1, 1)
arma$variance(10, 1)
arma$variance(30, 1)

## Update the coefficients
new_coef <- arma$coefficients * 1.01
arma$update(new_coef)
arma$coefficients == new_coef
sum(is.na(arma$std.errors)) == length(new_coef)

## Update the standard errors
new_std.errors <- arma$coefficients * 0.01
arma$update_std.errors(new_std.errors)
arma$std.errors == new_std.errors
sum(is.na(arma$std.errors)) == 0

## Update the variance sigma2
arma$sigma2
new_sigma2 <- 1.1
arma$update_sigma2(new_sigma2)
arma$sigma2 == new_sigma2

## Next step
x0 <- x[c(10)]
eps0 <- rnorm(1)
x_t <- list(`t+0` = c(x0, eps0))
# Simple forecast
x_t[[2]] <- arma$next_step(x_t[[1]], eps = 0)
# Check
x_t[[2]][1,1] == sum(arma$phi * x0) + sum(arma$theta * eps0)
# Two steps ahead
x_t[[3]] <- arma$next_step(x_t[[2]], eps = 0)
# Check
sum(arma$phi * x_t[[2]][1,1]) + sum(arma$theta * 0)  == x_t[[3]][1,1]






