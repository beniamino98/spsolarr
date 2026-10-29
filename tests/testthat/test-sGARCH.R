library(solarr)
library(rugarch)

x <- rnorm(1000)

# Initialize the model
GARCH <- sGARCH$new(1, 1)
# Fit
GARCH$fit(x)
# self <- GARCH$.__enclos_env__$self
# private <- GARCH$.__enclos_env__$private
GARCH

# Slots
# Order
GARCH$archOrder
GARCH$garchOrder
GARCH$order
# Parameters
GARCH$omega
GARCH$alpha
GARCH$beta
GARCH$coefficients
# Tidy output
GARCH$tidy
# Standard errors
GARCH$std.errors
# Matrices
GARCH$A
GARCH$b
GARCH$d
# Check unconditional vol
GARCH$sigma2_inf
# Total likelihood
GARCH$loglik
# Tidy output
GARCH$tidy

# Check next step prediction
eps0 <- c(0.3)
sigma20 <- c(0.9)
sigma2t <- sigma20
# One step ahead
next_step_garch <- function(omega = 1, alpha = 0, beta = 0, eps0, sigma20){
  GARCH$omega + GARCH$alpha * eps0^2 +GARCH$beta * sigma20
}
# One step ahead
sigma2t[2] <- next_step_garch(omega, alpha, beta, eps0, sigma2t[1])
sigma2t[2]
GARCH$next_step(eps0, sigma20)
# Two steps ahead
sigma2t[3] <- next_step_garch(omega, alpha, beta, sigma2t[2], sigma2t[2])
sigma2t[3]
GARCH$next_step(sigma2t[2], sigma2t[2])

# All sample prediction
GARCH$filter(x = x)


