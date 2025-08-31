rm(list = ls())

library(stats)
library(sinew)

# -----------------------------------------------------------------------------
# ---- Kernel Density Estimation with Mean Integrated Squared Error (MISE) ----
# -----------------------------------------------------------------------------

#' @title sigma_n
#' @description Gernerate the robust estimate of standard deviation
#' @param x A numeric input vector
#' @return a numeric value
#' @export
sigma_n <- function(x) {
  min(IQR(x) / (qnorm(0.75) - qnorm(0.25)), sd(x))
}

# Define functions for the normal pdf and its second derivative
#' @title phi
#' @description Create the normal pdf
#' @param u A numeric value
#' @return A density value at a given point u
#' @export
phi <- function(u) (1 / sqrt(2 * pi)) * exp(-0.5 * u^2)

#' @title phi_pp
#' @description Create the second derivative of the normal pdf
#' @param u A numeric value
#' @return A numeric value of the second derivative of the density function at a given point u
#' @export
phi_pp <- function(u) (u^2 - 1) * phi(u)

#' @title fpp_hat
#' @description Estimate f'' on a grid with the pilot bandwidth
#' @param x A numeric vector of x_i
#' @param y A numeric vector of x
#' @param h A numeric value of pilot bandwidth
#' @return The kernel-based estimate of curvature of the density function
#' @export
fpp_hat <- function(x, y, h) {
  u <- function(y, x) (y - x) / h
  u_mat <- outer(y, x, u)
  return(rowMeans(phi_pp(u_mat)) / h^3)
}

#' @title R_fpp_hat
#' @description Calculate the integrated squared second derivative of the density function.
#' @param x A numeric input vector
#' @param h A numeric value of bandwidth
#' @param multiplier A numeric multiplier for evaluation points
#' @param n_eval The number of KDE evaluation points
#' @return A numeric value of the integrated squared second derivative of pdf
#' @export
R_fpp_hat <- function(x, h, multiplier, n_eval) {
  lower_bound <- min(x) - multiplier * sigma_n(x)
  upper_bound <- max(x) + multiplier * sigma_n(x)
  y <- seq(lower_bound, upper_bound, length.out = n_eval)
  fpp_n <- fpp_hat(x, y, h)
  dx <- min(diff(y))
  return(sum(fpp_n^2 * dx))
}

#' @title update_bandwidth
#' @description Update the bandwidth of KDE
#' @param n The sample size
#' @param R_fpp_n The integrated squared second derivative
#' @return The updated bandwidth
#' @details DETAILS
#' @export
update_bandwidth <- function(n, R_fpp_n) {
  (1 / (2 * sqrt(pi) * R_fpp_n))^(1 / 5) * n^(-1 / 5)
}

#' @title bandwidth_mise
#' @description Find the optimal bandwidth
#' @param x A numeric input vector
#' @param guess_bandwidth The initial bandwidth
#' @param epsilon The difference between the old and updated bandwidths
#' @param multiplier A numeric multiplier for evaluation points
#' @param n_eval The number of KDE evaluation points
#' @return The converged bandwidth
#' @export
bandwidth_mise <- function(x, guess_bandwidth, epsilon, multiplier, n_eval) {
  n <- length(x)

  message(sprintf(
    "guess_bandwidth = %.5f", guess_bandwidth
  ))

  R_fpp_n <- R_fpp_hat(x, guess_bandwidth, multiplier, n_eval)
  bandwidth_update <- update_bandwidth(n, R_fpp_n)
  iteration <- 0

  while (abs(bandwidth_update - guess_bandwidth) > epsilon) {
    iteration <- iteration + 1
    guess_bandwidth <- bandwidth_update
    R_fpp_n <- R_fpp_hat(x, guess_bandwidth, multiplier, n_eval)
    bandwidth_update <- update_bandwidth(n, R_fpp_n)
    message(sprintf(
      "iteration %d  h_old = %.5f  bandwidth_update = %.5f",
      iteration, guess_bandwidth, bandwidth_update
    ))
  }

  return(bandwidth_update)
}

# Example
set.seed(13579)

epsilon <- 1e-6
multiplier <- 5
n_eval <- 1e3

x <- c(
  rnorm(1e4, mean = -4, sd = 0.8),
  rnorm(1e4, mean = -2, sd = 0.5),
  rnorm(1e4, mean = 4, sd = 1.0)
)

guess_bandwidth <- 5
bandwidth_s <- bandwidth_mise(x, guess_bandwidth, epsilon, multiplier, n_eval)
cat(sprintf("\nEstimated bandwidth: h = %.5f\n", bandwidth_s))

# KDE comparison
hist(x,
  breaks = 70,
  probability = TRUE, col = "lightgray",
  border = "white", main = "Kernel Density Estimation"
)
lines(density(x, bw = guess_bandwidth), col = "blue", lwd = 2)
lines(density(x, bw = bandwidth_s), col = "firebrick", lwd = 2)
legend("topright",
  c(
    sprintf("KDE (initial bandwidth): h=%.3f", guess_bandwidth),
    sprintf("KDE (optimal bandwidth): h=%.3f", bandwidth_s)
  ),
  col = c("blue", "firebrick"), lty = c(1, 2), lwd = 2, bty = "n"
)
