# VGAM prereleased dtruncpareto
dtruncpareto <-
  function(x, lower, upper, shape, log = FALSE) {
    # 20110324, log.arg added
    
    if (!is.logical(log.arg <- log) || length(log) != 1)
      stop("bad input for argument 'log'")
    rm(log)
    
    # if (!is.Numeric(x))
    #   stop("bad input for argument 'x'")
    if (!is.Numeric(lower, positive = TRUE))
      stop("argument 'lower' must be positive")
    if (!is.Numeric(upper, positive = TRUE))
      stop("argument 'upper' must be positive")
    if (!is.Numeric(shape, positive = TRUE))
      stop("argument 'shape' must be positive")
    
    L <- max(length(x), length(lower), length(upper), length(shape))
    if (length(x)     != L) x     <- rep_len(x,     L)
    if (length(shape) != L) shape <- rep_len(shape, L)
    if (length(lower) != L) lower <- rep_len(lower, L)
    if (length(upper) != L) upper <- rep_len(upper, L)
    
    logdensity <- rep_len(log(0), L)
    xok <- (0 <  lower) & (lower <= x) &
      (x <= upper) & (shape >  0)
    
    logdensity[xok] <- log(shape[xok]) +
      shape[xok] * log(lower[xok]) -
      (shape[xok] + 1) * log(x[xok]) -
      log1p(-(lower[xok] / upper[xok])^(shape[xok]))
    
    logdensity[shape <= 0] <- NaN
    logdensity[upper < lower] <- NaN
    logdensity[0 > lower] <- NaN
    
    if (log.arg) logdensity else exp(logdensity)
  }  # dtruncpareto


loglambertLambdaderiv <- function(x) {
  if (length(x) == 1) {
    output <- ifelse(x == 0, 0, log(x) - log(1 + x))
  } else {
    zeroindex <- which(x == 0)
    output <- log(x) - log(1 + x)
    output[zeroindex] <- 0
  }
}

###############
# joint prior
qgamma_transformed <- function(u1, u2, shape, rate) { # a function that returns a unique log(Gamma(shape, rate)) random variable based on two uniform random variables.
  # do NOT use this function for nested sampling - as the generation of a single gamma variate depends on two random variates, which doesn't obey pigeon hole and may confuse the algorithm
  # see https://math.stackexchange.com/questions/190670/how-exactly-are-the-beta-and-gamma-distributions-related
  # https://en.wikipedia.org/wiki/Beta_distribution#Related_distributions
  # https://en.wikipedia.org/wiki/Gamma_distribution
  # uniform to exponential: gexp ~ exp(shape)
  gexp <- qexp(u1, rate = shape)
  
  # exponential to beta: gbeta = exp(-gexp) ~ beta(shape, 1)
  #gbeta = exp(-gexp), so log(gbeta) = -gexp
  
  # beta to gamma: for V ~ gamma(1 + shape, rate), X = gbeta*V ~ gamma(shape, rate) and Y = (1-gbeta)*V ~ gamma(1, rate)
  lV <- log(qgamma(u2, 1 + shape, rate))
  Y <- lV + logminus(1, -gexp)
  X <- lV - gexp
  c(X, Y)
}

###############################
# broken power law
dbrokenpower <- function(x, pllocation_vec, plshape_vec) {
  pos_index <- ifelse(x >= max(pllocation_vec), # whether x is in truncated Pareto or Pareto
                      0, # x in Pareto => mark the upper limit as the 0th index
                      which(pllocation_vec - c(x) > 0)[1]) # truncated Pareto
  neg_index <- ifelse(pos_index != 0, pos_index - 1, length(pos_index))
  
  pllocation_j <- pllocation_vec[neg_index]
  plshape_j <- plshape_vec[neg_index]
  ifelse(pos_index == 0, VGAM::dpareto(x, scale = pllocation_j, shape = plshape_j, log = TRUE),
         dtruncpareto(x, lower = pllocation_j, upper = pllocation_vec[pos_index], shape = plshape_j, log = TRUE))
}
dbrokenpower(x = 3, pllocation_vec = 1:5, plshape_vec = 0.1 * (1:5))
x_brokenvec <- seq(1, 7, by = 0.1)
plot(x_brokenvec, sapply(x_brokenvec, function(i) dbrokenpower(x = i, pllocation_vec = 1:5, plshape_vec = 0.1 * (1:5))), type = "l")


#####################
# complete integrated kernel
########################
int_lik_test <- function(j) {
  D = last_py$simulated_data
  integrated_likelihood_transform(j[1], j[2], j[3], j[4], D, A, Time, a, R, e)
} 
#optim(c(mu, theta, pid, xi), fn = int_lik_test, control = list("fnscale" = -1), lower = rep(0, 4), upper = c(Inf, Inf, 1, Inf), method = "Nelder-Mead")
#library(dfoptim)
#nmkb(c(mu, theta, pid, xi), fn = int_lik_test, control = list("maximize" = TRUE), lower = rep(0, 4), upper = c(Inf, Inf, 1, Inf))

################
# transformed functions for dynesty
#################
int_prior_transform_emx <- function(u, x) { # must have some initial points that u[xi] > 0.999
  if (any(u < .Machine$double.xmin) || any(u >= 1)) {
    -Inf
  } else {
    muloc = 1.2 * 10^-6
    thetaloc = 8.7 * 10^-12
    lmu0 = log(10^6) - log(A) - log(Time)
    ltheta0 = log(10^18) - 2 * (log(A) + log(Time))
    
    umu <- log(LaplacesDemon::qtrunc(u[1], spec = "cauchy", location = muloc, scale = muloc * 100, a = 0))
    utheta <- log(LaplacesDemon::qtrunc(u[2], spec = "cauchy", location = thetaloc, scale = thetaloc * 100, a = 0))
    upid <- LaplacesDemon::logit(u[3])
    uxi <- log(qgamma(u[4], shape = x, rate = 1))
    #utxi <- qgamma_transformed(u[4], u[5], shape = exp(2 * lmu0 - ltheta0), rate = exp(lmu0 - ltheta0))
    
    c(umu, utheta, upid, uxi)
  }
}

prior_transform_emx <- function(u, x) { # must have some initial points that u[xi] > 0.999
  if (any(u < .Machine$double.xmin) || any(u >= 1)) {
    -Inf
  } else {
    muloc = 1.2 * 10^-6
    thetaloc = 8.7 * 10^-12
    lmu0 = log(10^6) - log(A) - log(Time)
    ltheta0 = log(10^18) - 2 * (log(A) + log(Time))
    
    umu <- log(LaplacesDemon::qtrunc(u[1], spec = "cauchy", location = muloc, scale = muloc * 100, a = 0))
    utheta <- log(LaplacesDemon::qtrunc(u[2], spec = "cauchy", location = thetaloc, scale = thetaloc * 100, a = 0))
    upid <- LaplacesDemon::logit(u[3])
    uxi <- log(qgamma(u[4], shape = x, rate = 1))
    #utxi <- qgamma_transformed(u[4], u[5], shape = exp(2 * lmu0 - ltheta0), rate = exp(lmu0 - ltheta0))
    
    uLambda <- q0inf_gamma(u[-(1:4)], shape = exp(2 * umu - utheta), rate = exp(umu - utheta), pid = u[3])
    
    c(umu, utheta, upid, uxi, 0, uLambda)
  }
}

##################
# posterior kernel
##################
logpostker_emx <- function(mu, theta, pid, xi, Lambda, 
                           D, A, Time, a, r, e,
                           muloc = 1.2 * 10^-6, thetaloc = 8.7 * 10^-12)
  loglik(xi, Lambda, D, A, Time, a, r, e) +
  logxiprior_xi(xi, D[length(D)]) +
  logLambdaprior_gamma(Lambda, mu, theta, pid) +
  logpidhypprior(pid) +
  logmuhypprior(mu, muloc = muloc) +
  logthetahypprior(theta, thetaloc = thetaloc)

t_logpostker_emx <- function(t_mu, t_theta, t_pid, t_xi, t_Lambda, 
                             D, A, Time, a, r, e, 
                             muloc = 1.2 * 10^-6, thetaloc = 8.7 * 10^-12) {
  t_loglik(t_xi, t_Lambda, D, A, Time, a, r, e) +
    t_logxiprior_gamma(t_xi, D[length(D)]) +
    t_logLambdaprior_gamma(t_Lambda, t_mu, t_theta, t_pid) +
    t_logpidhypprior(t_pid) +
    t_logmuhypprior(t_mu, muloc = muloc) +
    t_logthetahypprior(t_theta, thetaloc = thetaloc)
}
