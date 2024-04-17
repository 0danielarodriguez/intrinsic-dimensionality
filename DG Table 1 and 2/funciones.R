elije_donoho_el_d<-function(Y){
  
# Assuming you have a matrix Y
# And assuming you have the function optimal_SVHT_coef_sigma_unknown available in R

# Perform Singular Value Decomposition
svd_result <- svd(scale(Y,center=TRUE,scale=FALSE))

# Extract U, D, and V matrices
U <- svd_result$u
D <- svd_result$d
V <- svd_result$v

# Diagonal values of Y, but in MATLAB 'diag(Y)' might be a mistake if Y is meant.
# Assuming it is meant to be diag(D) which are the singular values of Y.
y <- D

# Apply thresholding
threshold <- optimal_SVHT_coef_sigma_unknown(min(c(nrow(Y)/ncol(Y), ncol(Y)/nrow(Y)))) * median(y) 
y[y < threshold] <- 0

# Construct denoised matrix Xhat
Xhat <- U %*% diag(y) %*% t(V)

# Xhat is the denoised version of Y
d=sum(y>threshold)
d
}



optimal_SVHT_coef_sigma_known <- function(beta) {
  stopifnot(all(beta > 0))
  stopifnot(all(beta <= 1))
  stopifnot(is.vector(beta))
  
  w <- (8 * beta) / (beta + 1 + sqrt(beta^2 + 14 * beta + 1))
  lambda_star <- sqrt(2 * (beta + 1) + w)
  
  return(lambda_star)
}


optimal_SVHT_coef_sigma_unknown <- function(beta) {
  # No direct equivalent to MATLAB's warning off in R, but generally, warnings can be suppressed or managed globally
  stopifnot(all(beta > 0))
  stopifnot(all(beta <= 1))
  stopifnot(length(beta) == length(beta)) # Ensuring beta is a vector
  
  coef <- optimal_SVHT_coef_sigma_known(beta)
  
  MPmedian <- numeric(length(beta))
  for (i in seq_along(beta)) {
    MPmedian[i] <- MedianMarcenkoPastur(beta[i])
  }
  
  omega <- coef / sqrt(MPmedian)
  return(omega)
}

MarcenkoPasturIntegral <- function(x, beta) {
  if (beta <= 0 | beta > 1) {
    stop('beta beyond')
  }
  lobnd <- (1 - sqrt(beta))^2
  hibnd <- (1 + sqrt(beta))^2
  if (x < lobnd | x > hibnd) {
    stop('x beyond')
  }
  dens <- function(t) sqrt((hibnd-t)*(t-lobnd)) / (2*pi*beta*t)
  I <- integrate(dens, lobnd, x)$value
  cat(sprintf('x=%.3f,beta=%.3f,I=%.3f\n', x, beta, I))
  return(I)
}


MedianMarcenkoPastur <- function(beta) {
  MarPas <- function(x) 1 - incMarPas(x, beta, 0)
  lobnd <- (1 - sqrt(beta))^2
  hibnd <- (1 + sqrt(beta))^2
  change <- TRUE
  while (change & (hibnd - lobnd > .001)) {
    change <- FALSE
    x <- seq(lobnd, hibnd, length.out = 5)
    y <- numeric(length(x))
    for (i in seq_along(x)) {
      y[i] <- MarPas(x[i])
    }
    if (any(y < 0.5)) {
      lobnd <- max(x[y < 0.5])
      change <- TRUE
    }
    if (any(y > 0.5)) {
      hibnd <- min(x[y > 0.5])
      change <- TRUE
    }
  }
  med <- (hibnd + lobnd) / 2
  return(med)
}


incMarPas <- function(x0, beta, gamma) {
  if (beta > 1) {
    stop('betaBeyond')
  }
  topSpec <- (1 + sqrt(beta))^2
  botSpec <- (1 - sqrt(beta))^2
  MarPas <- function(x) ifelse((topSpec-x)*(x-botSpec) > 0,
                               sqrt((topSpec-x)*(x-botSpec)) / (beta * x) / (2 * pi),
                               0)
  if (gamma != 0) {
    fun <- function(x) (x^gamma * MarPas(x))
  } else {
    fun <- MarPas
  }
  I <- integrate(fun, x0, topSpec)$value
  return(I)
}


incMarPas <- function(x0, beta, gamma) {
  if (beta > 1) {
    stop('betaBeyond')
  }
  topSpec <- (1 + sqrt(beta))^2
  botSpec <- (1 - sqrt(beta))^2
  MarPas <- function(x) ifelse((topSpec-x)*(x-botSpec) > 0,
                               sqrt((topSpec-x)*(x-botSpec)) / (beta * x) / (2 * pi),
                               0)
  if (gamma != 0) {
    fun <- function(x) (x^gamma * MarPas(x))
  } else {
    fun <- MarPas
  }
  I <- integrate(fun, x0, topSpec)$value
  return(I)
}



