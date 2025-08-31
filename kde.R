rm(list = ls())

library(stats)

# Kernel Density Estimation with Mean Integrated Squared Error (MISE)
sigma_n <- function(x) {
  min(IQR(x) / (qnorm(0.75) - qnorm(0.25)), sd(x))
}

# Initialize bandwidth
initial_bandwidth <- function(x, rho) rho * sigma_n(x) * length(x)^(-1 / 5)

# Define functions for the normal pdf and its second derivative
phi <- function(u) (1 / sqrt(2 * pi)) * exp(-0.5 * u^2)
phi_pp <- function(u) (u^2 - 1) * phi(u)

# Estimate f''
fpp_hat <- function(x, y, h, k) {
  u <- function(y, x) (y - x) / h
  U <- outer(y, x, u)
  return(rowMeans(phi_pp(U)) / (h^k))
}

# Estimate R(f'')
R_fpp_hat <- function(x, h, multiplier, m) {
  lower_bound <- min(x) - multiplier * sigma_n(x)
  upper_bound <- max(x) + multiplier * sigma_n(x)
  y <- seq(lower_bound, upper_bound, length.out = m)
  fpp_n <- fpp_hat(x, y, h, k = 1)
  dx <- min(diff(y))
  return(sum(fpp_n^2 * dx))
}

# Update bandwidth
bandwidth_update <- function(n, R_fpp_n) {
  ((1 / (2 * sqrt(pi)) / R_fpp_n)^(1 / 5)) * n^(-1 / 5)
}

# Find the optimal bandwidth
bandwidth_mise <- function(x, epsilon, multiplier, m) {
  n <- length(x)
  h_0 <- initial_bandwidth(x, rho = 0.9)

  R_fpp_n <- R_fpp_hat(x, h_0, multiplier, m)
  h_update <- bandwidth_update(n, R_fpp_n)
  iteration <- 0

  while (abs(h_update - h_0) > epsilon) {
    iteration <- iteration + 1
    h_0 <- h_update
    R_fpp_n <- R_fpp_hat(x, h_0, multiplier, m)
    h_update <- bandwidth_update(n, R_fpp_n)
    message(sprintf(
      "iteration=%d  h_old=%.5f  h_update=%.5f",
      iteration, h_0, h_update
    ))
  }

  return(h_0)
}

# Example
set.seed(13579)
x <- c(
  rnorm(1e4, mean = -4, sd = 0.8),
  rnorm(1e4, mean = -2, sd = 0.5),
  rnorm(1e4, mean = 4, sd = 1.0)
)

h_0 <- initial_bandwidth(x, rho = 5)
h_s <- bandwidth_mise(x, epsilon = 1e-6, multiplier = 5, m = 5e3)
cat(sprintf("\nEstimated bandwidth: h = %.5f\n", h_s))

# KDEs comparison
hist(x,
  breaks = 70,
  probability = TRUE, col = "lightgray",
  border = "white", main = "Kernel Density Estimation",
  xlab = "x"
)
lines(density(x, bw = h_0), col = "blue", lwd = 2)
lines(density(x, bw = h_s), col = "firebrick", lwd = 2)
legend("topright",
  c(
    sprintf("KDE (initial bandwidth): h=%.3f", h_0),
    sprintf("KDE (optimal bandwidth): h=%.3f", h_s)
  ),
  col = c("blue", "firebrick"), lty = c(1, 2), lwd = 2, bty = "n"
)
