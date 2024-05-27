
# Function for random sampling of the generalized inverse Gaussian distribution
# Implementation from Devroye (2014) algorithm
gigrnd <- function(lambda, a, b, sample_size) {
  
  # Setup: sample from 2 parameter version of GIG(alpha, omega)
  omega <- sqrt(a * b)
  
  swap <- FALSE
  if (lambda < 0) {
    lambda <- -lambda
    swap <- TRUE
  }
  
  alpha <- sqrt(omega^2 + lambda^2) - lambda
  
  # Find t
  x <- -psi(1, alpha, lambda)
  if ((x >= 0.5) && (x <= 2)) {
    t <- 1
  }
  else if (x > 2) {
    t <- sqrt(2 / (alpha + lambda))
  }
  else {
    t <- log(4 / (alpha + 2 * lambda))
  }
  
  # Find s
  x <- -psi(1, alpha, lambda)
  if ((x >= 0.5) && (x <= 2)) {
    s <- 1
  }
  else if (x > 2) {
    s <- sqrt(4 / (alpha * cosh(1) + lambda))
  }
  else {
    s <- log(4 / (alpha + 2 * lambda))
  }
  
  # Generation
  eta <- -psi(t, alpha, lambda)
  zeta <- -dpsi(t, alpha, lambda)
  theta <- -psi(-s, alpha, lambda)
  xi <- dpsi(-s, alpha, lambda)
  p <- 1 / xi
  r <- 1 / zeta
  td <- t - (r * eta)
  sd <- s - (p * theta)
  q <- td + sd

  X <- numeric(sample_size)
  for (sample in 1:sample_size) {
    done <- FALSE
    while (!done) {
      U <- runif(1)
      V <- runif(1)
      W <- runif(1)
      if (U < (q / (p + q + r))) {
        X[sample] <- -sd + (q * V)
      }
      else if (U < ((q + r) / (p + q + r))) {
        X[sample] <- td - (r * log(V))
      }
      else {
        X[sample] <- -sd + (p * log(V))
      }
      
      # Check to terminate
      f1 <- exp(-eta - zeta * (X[sample] - t))
      f2 <- exp(-theta + xi * (X[sample] + s))
      if ((W * fun_g(X[sample], sd, td, f1, f2)) <= 
          exp(psi(X[sample], alpha, lambda))) {
        done <- TRUE
      }
    }
  }
  
  X <- exp(X) * (lambda / omega + sqrt(1 + (lambda / omega)^2))
  if (swap) {
    X <- 1 / X
  }
  
  return(X / sqrt(a / b))
}

psi <- function(x, alpha, lambda) {
  return(-alpha * (cosh(x) - 1) - lambda * (exp(x) - x - 1))
}

dpsi <- function(x, alpha, lambda) {
  return(-alpha * sinh(x) - lambda * (exp(x) - 1))
}

fun_g <- function(x, sd, td, f1, f2) {

  rval <- 0
  if ((x >= -sd) && (x <= td)) {
    rval <- 1
  }
  else if (x > td) {
    rval <- f1
  }
  else if (x < -sd) {
    rval <- f2
  }
  
  return(rval)
}









gigrnd_debug <- function(lambda, a, b, sample_size) {
  
  # Setup: sample from 2 parameter version of GIG(alpha, omega)
  omega <- sqrt(a * b)
  
  swap <- FALSE
  if (lambda < 0) {
    lambda <- -lambda
    swap <- TRUE
  }
  
  alpha <- sqrt(omega^2 + lambda^2) - lambda
  
  # Find t
  x <- -psi(1, alpha, lambda)
  if ((x >= 0.5) && (x <= 2)) {
    t <- 1
  }
  else if (x > 2) {
    t <- sqrt(2 / (alpha + lambda))
  }
  else {
    t <- log(4 / (alpha + 2 * lambda))
  }
  
  # Find s
  x <- -psi(1, alpha, lambda)
  if ((x >= 0.5) && (x <= 2)) {
    s <- 1
  }
  else if (x > 2) {
    s <- sqrt(4 / (alpha * cosh(1) + lambda))
  }
  else {
    s <- log(4 / (alpha + 2 * lambda))
  }
  
  # Generation
  eta <- -psi(t, alpha, lambda)
  zeta <- -dpsi(t, alpha, lambda)
  theta <- -psi(-s, alpha, lambda)
  xi <- dpsi(-s, alpha, lambda)
  p <- 1 / xi
  r <- 1 / zeta
  td <- t - (r * eta)
  sd <- s - (p * theta)
  q <- td + sd
  
  usample <- c()
  X <- numeric(sample_size)
  for (sample in 1:sample_size) {
    done <- FALSE
    while (!done) {
      U <- runif(1)
      V <- runif(1)
      W <- runif(1)
      cat("R CODE U V W: ", U, V, W, "\n")
      usample <- c(usample, U, V, W)
      if (U < (q / (p + q + r))) {
        cat("case1\n")
        X[sample] <- -sd + (q * V)
      }
      else if (U < ((q + r) / (p + q + r))) {
        cat("case2\n")
        X[sample] <- td - (r * log(V))
      }
      else {
        cat("case3\n")
        X[sample] <- -sd + (p * log(V))
      }
      cat("sample generated in loop", X[sample], "\n")
      # Check to terminate
      f1 <- exp(-eta - zeta * (X[sample] - t))
      f2 <- exp(-theta + xi * (X[sample] + s))
      if ((W * fun_g(X[sample], sd, td, f1, f2)) <= 
          exp(psi(X[sample], alpha, lambda))) {
        done <- TRUE
      }
    }
  }
  
  X <- exp(X) * (lambda / omega + sqrt(1 + (lambda / omega)^2))
  if (swap) {
    cat("did swap\n")
    X <- 1 / X
  }
  
  return(list(ans=X / sqrt(a / b), usample=usample))
}