library(zipfR)
library(Rcpp)
library(VGAM) # (truncated) broken power law
library(LaplacesDemon) # for logistic function
library(lamW)
sourceCpp("logsum.cpp")

#############
# likelihood
#############
loglik <- function(xi, Lambda, 
                   D, A, Time, a, r, e) {
  #if (length(Lambda) != length(D) - 1){
  #  print(paste("Lambda:", Lambda))
  #  print(paste("D:", D))
  #  stop("Lambda and D do not have the same length.")
  #}
  if (xi >= .Machine$double.xmin) {
    dpois(x = D[length(D)],
          A * Time * xi, 
          log = TRUE) +
      sum(dpois(x = D[-length(D)], 
                (c(a * xi) + c(r) * c(e) * Lambda) * c(Time), 
                log = TRUE))
  } else { # large negative finite number for dynesty
    as.numeric(-1e100)
  }
}


t_loglik <- function(t_xi, t_Lambda, D, A, Time, a, r, e) {
  xi <- exp(t_xi)
  #Lambda <- lambertW0(exp(t_Lambda))
  Lambda <- t_Lambda
  
  dpois(x = D[length(D)],
        A * Time * xi, 
        log = TRUE) +
    sum(dpois(x = D[-length(D)], 
              (c(a * xi) + c(r) * c(e) * Lambda) * c(Time), 
              log = TRUE))
}

################
# integrated likelihood
################
source(paste0(getwd(), "/numerical_integration.R"))
integrated_kernel <- function(t_mu = NULL, t_theta = NULL, 
                              t_shape = 2 * t_mu - t_theta, 
                              t_rate = t_mu - t_theta,
                              t_pid, t_xi, 
                              D, n = length(D) - 1, A, Time, a, r, e,
                              shape1 = 1, shape2 = 1, muloc = 1.2e-6, thetaloc = 8.7e-12, 
                              mu0 = 10^6 / (A * Time), theta0 = 10^18 / (A * Time)^2,
                              analytical_method = TRUE, Lamba_seq = seq(0, exp(7), by = 0.1)) {
  # integral function: using analytical or numerical:
  integral <- ifelse(analytical_method,
                     analytical(t_shape = t_shape, t_rate = t_rate, t_xi = t_xi, t_pid = t_pid, D = D, A = A, Time = Time, a = a, r = r, e = e),
                     numerical(t_mu, t_theta, t_pid, t_xi,
                               Lambda = Lambda_seq,
                               x = D[n + 1], y = D[1:n], A, Time, a, R, e = e))
  
  # other functions
  hyperpriors <- t_logpidhypprior(t_pid, shape1, shape2) +
    t_logshaperatehypprior(t_shape, t_rate, muloc, thetaloc) +
    t_logxiprior_gamma(t_xi, mu0, theta0)
  if (!is.null(t_mu) && !is.null(t_theta)) {
    hyperpriors <- t_logpidhypprior(t_pid, shape1, shape2) +
      t_logmuhypprior(t_mu) +
      t_logthetahypprior(t_theta) +
      t_logxiprior_gamma(t_xi, mu0, theta0)
  }

  integral + hyperpriors
}


t_int_kernel_gamma <- function(obj, D, A, Time, a, R, e, mu_theta_param = TRUE)
  UseMethod("t_int_kernel_gamma")

t_int_kernel_gamma.numeric <- function(zeta, D, A, Time, a, R, e, mu_theta_param = TRUE) {
  ll <- ifelse(mu_theta_param,
               integrated_kernel(t_mu = zeta[1],
                                 t_theta = zeta[2],
                                 t_pid = zeta[3],
                                 t_xi = zeta[4],
                                 D = D,
                                 A = A,
                                 Time = Time,
                                 a = a,
                                 r = R,
                                 e = e),
               integrated_kernel(t_shape = zeta[1],
                                 t_rate = zeta[2],
                                 t_pid = zeta[3],
                                 t_xi = zeta[4],
                                 D = D,
                                 A = A,
                                 Time = Time,
                                 a = a,
                                 r = R,
                                 e = e))
  
  ifelse(is.na(ll) || is.infinite(ll), -Inf, ll) # in case the function returns NaN, reject that point
}

t_int_kernel_gamma.matrix <- function(zeta, D, A, Time, a, R, e, mu_theta_param = TRUE) {
  apply(zeta, 1, t_int_kernel_gamma.numeric,
        D = D, A = A, Time = Time, a = a, R = R, e = e, mu_theta_param = mu_theta_param)
}

#######
# prior
#######
#####################
# gamma
logdeltafunction <- function(i) ifelse(i == 0, 0, -Inf)

logLambdaprior_gamma <- function(Lambda,
                                 mu, theta, pid) {
  sapply(Lambda, 
         function(i) logplus(log(pid) + logdeltafunction(i),
                             log(1 - pid) + dgamma(x = i,
                                                   shape = mu^2 / theta,
                                                   rate = mu / theta,
                                                   log = TRUE))) |> sum()
}
  
logxiprior_gamma <- function(xi,
                             mu0, theta0)
  dgamma(xi, 
         shape = mu0^2 / theta0,
         rate = mu0 / theta0,
         log = TRUE)

logxiprior_emx <- function(xi, x)
  dgamma(xi, 
         shape = x,
         rate = 1,
         log = TRUE)

dgamma_infcorr <- function(x, shape, rate, log) {
  ll <- dgamma(x, shape, rate = rate, log = log)
  if (is.infinite(ll)) {
    return(-Inf)
  } else {
    return(ll)
  }
}

t_logLambdaprior_gamma <- function(t_Lambda, # Lambda+log(Lambda)
                                   t_mu, # log(mu)
                                   t_theta, # log(theta)
                                   t_pid # LaplacesDemon::logit(pid)
) {
  #Lambda <- lambertW0(exp(t_Lambda))
  Lambda <- t_Lambda
  mu <- exp(t_mu)
  theta <- exp(t_theta)
  pid <- invlogit(t_pid)
  # log Jacobian for : sum(loglambertLambdaderiv(x = Lambda))
  sapply(Lambda, # Jacobian: exp(sum(t_Lambda))
                         function(i) logplus(log(pid) + logdeltafunction(i),
                                             log(1 - pid) + dgamma_infcorr(x = i,
                                                                   shape = exp(2 * t_mu - t_theta),
                                                                   rate = exp(t_mu - t_theta),
                                                                   log = TRUE))) |> sum() #Jacobian:prod(Lambda) => log Jacobian:sum(t_Lambda)
}

t_logLambdaprior_puregamma <- function(t_Lambda, # Lambda+log(Lambda)
                                   t_mu, # log(mu)
                                   t_theta, # log(theta)
                                   t_pid # LaplacesDemon::logit(pid)
) {
  #Lambda <- lambertW0(exp(t_Lambda))
  Lambda <- t_Lambda
  mu <- exp(t_mu)
  theta <- exp(t_theta)
  pid <- invlogit(t_pid)
  sum(dgamma(x = Lambda, shape = exp(2 * t_mu - t_theta),
                     rate = exp(t_mu - t_theta), log = TRUE))
}

t_logxiprior_gamma <- function(t_xi, # log(xi)
                               mu0, theta0) {
  xi <- exp(t_xi)
  dgamma(xi, 
         shape = mu0^2 / theta0,
         rate = mu0 / theta0,
         log = TRUE) + t_xi # Jacobian: xi => log Jacobian: t_xi
}
t_logxiprior_emx <- function(t_xi, # log(xi)
                             x) {
  xi <- exp(t_xi)
  dgamma(xi, 
         shape = x,
         rate = 1,
         log = TRUE) + t_xi # Jacobian: xi => log Jacobian: t_xi
}



############
# hyperprior
############
logpidhypprior <- function(pid, shape1 = 1, shape2 = 1)
  dbeta(pid, shape1 = shape1, shape2 = shape2, log = TRUE)

logmuhypprior <- function(mu, muloc = 1.2e-6)
  ifelse(mu < 0,
         -Inf,
         dcauchy(mu, location = muloc, scale = 100 * muloc, log = TRUE))

logthetahypprior <- function(theta, thetaloc = 8.7e-12)
  ifelse(theta < 0,
         -Inf,
         dcauchy(theta, location = thetaloc, scale = 100 * thetaloc, log = TRUE))

t_logpidhypprior <- function(t_pid, # LaplacesDemon::logit(pid)
                             shape1 = 1, shape2 = 1) {
  pid <- invlogit(t_pid)
  dbeta(pid, shape1 = shape1, shape2 = shape2, log = TRUE) + log(pid) + logplus(0, log(pid) - t_pid)
}


t_logmuhypprior <- function(t_mu, muloc = 1.2 * 10^-6) {
  mu <- exp(t_mu)
  dcauchy(mu, location = muloc, scale = 100 * muloc, log = TRUE) + t_mu #Jacobian:mu; log Jacobian:t_mu
}

t_logthetahypprior <- function(t_theta, thetaloc = 8.7 * 10^-12) {
  theta <- exp(t_theta)
  dcauchy(theta, location = thetaloc, scale = 100 * thetaloc, log = TRUE) + t_theta #Jacobain:theta; log Jacobian:t_theta
}

t_logshaperatehypprior <- function(t_shape, t_rate, muloc = 1.2 * 10^-15, thetaloc = 8.7 * 10^-12) {
  mu <- exp(t_shape - t_rate)
  theta <- exp(t_shape - 2 * t_rate)
  
  dcauchy(mu, location = muloc, scale = 100 * muloc, log = TRUE) + 
    dcauchy(theta, location = thetaloc, scale = 100 * thetaloc, log = TRUE) # Jacobian is 1.
  # easy to interpret why Jacobian=1: the mass of the distribution distributes the same way on R^2, regardless of (shape, rate) or (mu, theta).
}
t_logshaperatehypprior <- function(t_shape, t_rate, muloc = 1.2 * 10^-6, thetaloc = 8.7 * 10^-12) {
  shapeloc <- muloc^2 / thetaloc
  scaleloc <- muloc / thetaloc
  dcauchy(exp(t_shape), location = shapeloc, scale = 100 * shapeloc, log = TRUE) + 
    dcauchy(exp(t_rate), location = scaleloc, scale = 100 * scaleloc, log = TRUE) + # Jacobian
    t_shape + t_rate
}

#############
# joint prior
#############
# see https://statproofbook.github.io/P/gam-qf.html and https://stats.stackexchange.com/questions/145262/zero-inflated-gamma-how-to-write-down-the-cdf
q0inf_gamma <- function(prob_vec, shape, rate, pid) {
  output_vec <- rep(NA, length(prob_vec))
  output_vec[prob_vec <= pid] <- 0
  
  gamma_part_p <- c(prob_vec[prob_vec > pid]) - pid
  trans_p <- lgamma(shape) + log(gamma_part_p) - log(1 - pid)
  output <- log(Igamma.inv(shape, trans_p, log = TRUE)) - log(rate)
  output_vec[prob_vec > pid] <- exp(output)
  
  output_vec
}
qgamma_bounded <- function(p, shape, rate) {
  plowbound <- pgamma(q = .Machine$double.eps, shape = shape, rate = rate)
  # transform the input p so it is bounded below by machine xmin
  pt <- (1 - plowbound) * p + plowbound
  qgamma(pt, shape, rate)
}
rgamma_transformed <- function(n, shape, rate) { # a function to sample from gamma distributions with tiny shape and rate on log scale.
  # do NOT use this function for nested sampling - as the generation of a single gamma variate depends on two random variates, which doesn't obey pigeon hole and may confuse the algorithm
  # see https://math.stackexchange.com/questions/190670/how-exactly-are-the-beta-and-gamma-distributions-related
  # https://en.wikipedia.org/wiki/Beta_distribution#Related_distributions
  # https://en.wikipedia.org/wiki/Gamma_distribution
  # uniform to exponential: gexp ~ exp(shape)
  gexp <- qexp(runif(n), rate = shape)
  
  # exponential to beta: gbeta = exp(-gexp) ~ beta(shape, 1)
  #gbeta = exp(-gexp), so log(gbeta) = -gexp
  
  # beta to gamma: for V ~ gamma(1 + shape, rate), X = gbeta*V ~ gamma(shape, rate) and Y = (1-gbeta)*V ~ gamma(1, rate)
  lV <- log(rgamma(n, 1 + shape, rate))
  Y <- lV + sapply(-gexp, function(j) logminus(1, j))
  X <- lV - gexp
  c(X, Y)
}


################
# transformed functions for dynesty
#################
# prior
int_prior_transform <- function(u) { # must have some initial points that u[xi] > 0.999
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
    
    uxi <- log(qgamma_bounded(u[4], shape = exp(2 * lmu0 - ltheta0), rate = exp(lmu0 - ltheta0)))
    #utxi <- qgamma_transformed(u[4], u[5], shape = exp(2 * lmu0 - ltheta0), rate = exp(lmu0 - ltheta0))
    
    c(umu, utheta, upid, uxi)
  }
  
}

prior_transform <- function(u) { # must have some initial points that u[xi] > 0.999
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
    uxi <- log(qgamma_bounded(u[4], shape = exp(2 * lmu0 - ltheta0), rate = exp(lmu0 - ltheta0)))
    #utxi <- qgamma_transformed(u[4], u[5], shape = exp(2 * lmu0 - ltheta0), rate = exp(lmu0 - ltheta0))
    
    #uLambda <- qgamma(u[-(1:4)], shape = exp(2 * log(umu) - log(utheta)), rate = exp(log(umu) - log(utheta)))
    uLambda <- q0inf_gamma(u[-(1:4)], shape = exp(2 * umu - utheta), rate = exp(umu - utheta), pid = u[3])
    
    c(umu, utheta, upid, uxi, uLambda)
  }
}

# NEW!
prior_over_transform <- function(u, K) { # must have some initial points that u[xi] > 0.999
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
    uxi <- log(qgamma_bounded(u[4:(3 + K)], shape = exp(2 * lmu0 - ltheta0), rate = exp(lmu0 - ltheta0)))
    #utxi <- qgamma_transformed(u[4], u[5], shape = exp(2 * lmu0 - ltheta0), rate = exp(lmu0 - ltheta0))
    
    #uLambda <- qgamma(u[-(1:4)], shape = exp(2 * log(umu) - log(utheta)), rate = exp(log(umu) - log(utheta)))
    uLambda <- q0inf_gamma(u[-(1:(3 + K))], shape = exp(2 * umu - utheta), rate = exp(umu - utheta), pid = u[3])
    
    c(umu, utheta, upid, uxi, uLambda)
  }
}

int_prior_over_transform <- function(u) { # must have some initial points that u[xi] > 0.999
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
    
    uxi <- log(qgamma_bounded(u[-(1:3)], shape = exp(2 * lmu0 - ltheta0), rate = exp(lmu0 - ltheta0)))
    #utxi <- qgamma_transformed(u[4], u[5], shape = exp(2 * lmu0 - ltheta0), rate = exp(lmu0 - ltheta0))
    
    c(umu, utheta, upid, uxi)
  }
  
}

# likelihood
log_dpois <- function(loglambda, x) {
  x * c(loglambda) - lfactorial(x) - c(exp(loglambda))
}
loglik_transform <- function(t_xi, Lambda, 
                   D, A, Time, a, r, e) {
  if (t_xi >= -.Machine$double.xmax && t_xi <= .Machine$double.xmax) {
    lfreq <- log(A) + log(Time) + t_xi
    lfreqa <- c(Time) * sapply(c(log(r)) + c(log(e)) + log(Lambda), logplus, y = log(a) + t_xi)
    log_dpois(loglambda = lfreq, 
              x = D[length(D)]) +
      sum(log_dpois(loglambda = lfreqa,
                    x = D[-length(D)]))# +
      #dnorm(x = aux, mean = aux, log = TRUE)
  } else { # large negative finite number for dynesty
    as.numeric(-1e100)
  }
}

integrated_likelihood_transform <- function(t_mu, t_theta, t_pid, t_xi, D,
                                  A, Time, a, R, e) {
  
  if (t_xi <= -.Machine$double.xmax || t_xi >= .Machine$double.xmax) {
    as.numeric(-1e100)
  } else {
    n = length(D) - 1
    shape1 = 1
    shape2 = 1
    muloc = 1.2e-6
    thetaloc = 8.7e-12
    mu0 = 10^6 / (A * Time)
    theta0 = 10^18 / (A * Time)^2
    #Lamba_seq = seq(0, exp(7), by = 0.1)
    
    analytical(t_mu = t_mu, t_theta = t_theta, t_xi = t_xi, t_pid = t_pid, D = D, A = A, Time = Time, a = a, r = R, e = e)
  }
}


##################
# posterior kernel
##################
logpostker_gamma <- function(mu, theta, pid, xi, Lambda, 
                             D, A, Time, a, r, e, 
                             mu0 = 10^6 / (A * Time), theta0 = 10^18 / (A * Time)^2,
                             muloc = 1.2 * 10^-6, thetaloc = 8.7 * 10^-12)
  loglik(xi, Lambda, D, A, Time, a, r, e) +
  logxiprior_gamma(xi, mu0, theta0) +
  logLambdaprior_gamma(Lambda, mu, theta, pid) +
  logpidhypprior(pid) +
  logmuhypprior(mu, muloc = muloc) +
  logthetahypprior(theta, thetaloc = thetaloc)

t_logpostker_gamma <- function(t_mu, t_theta, t_pid, t_xi, t_Lambda, 
                               D, A, Time, a, r, e, 
                               mu0 = 10^6 / (A * Time), theta0 = 10^18 / (A * Time)^2,
                               muloc = 1.2 * 10^-6, thetaloc = 8.7 * 10^-12) {
  t_loglik(t_xi, t_Lambda, D, A, Time, a, r, e) +
    t_logxiprior_gamma(t_xi, mu0, theta0) +
    t_logLambdaprior(t_Lambda, t_mu, t_theta, t_pid) +
    t_logpidhypprior(t_pid) +
    t_logmuhypprior(t_mu, muloc = muloc) +
    t_logthetahypprior(t_theta, thetaloc = thetaloc)
}

# S3 construction
logpostker_unified_gamma <- suppressWarnings(function(obj, D, A, Time, a, R, e)
  UseMethod("logpostker_unified_gamma"))

logpostker_unified_gamma.numeric <- function(zeta, D, A, Time, a, R, e)  {
  l = logpostker_gamma(mu = zeta[1],
                       theta = zeta[2],
                       pid = zeta[3],
                       xi = zeta[4],
                       Lambda = zeta[-(1:4)],
                       D = D,
                       A = A,
                       Time = Time,
                       a = a,
                       r = R, 
                       e = e)
  ifelse(is.na(l), -Inf, l)
}

logpostker_unified_gamma.matrix <- function(zeta, D, A, Time, a, R, e) {
  apply(zeta, 1, logpostker_unified_gamma.numeric,
        D = D, A = A, Time = Time, a = a, R = R, e = e)
}

t_logpostker_unified_gamma <- suppressWarnings(function(obj, D, A, Time, a, R, e)
  UseMethod("t_logpostker_unified_gamma"))

t_logpostker_unified_gamma.numeric <- function(zeta, D, A, Time, a, R, e) {
  l <- t_logpostker_gamma(t_mu = zeta[1],
                          t_theta = zeta[2],
                          t_pid = zeta[3],
                          t_xi = zeta[4],
                          t_Lambda = zeta[-(1:4)],
                          D = D,
                          A = A,
                          Time = Time,
                          a = a,
                          r = R, 
                          e = e)
  ifelse(is.na(l) || is.infinite(l), -Inf, l) # in case the function returns NaN, reject that point
}

t_logpostker_unified_gamma.matrix <- function(zeta, D, A, Time, a, R, e) {
  apply(zeta, 1, t_logpostker_unified_gamma.numeric,
        D = D, A = A, Time = Time, a = a, R = R, e = e)
}

################
# NEW! likelihood extension
#########
loglik_over <- function(xi, # vector of length K
                       Lambda, # vector of length I (i indexed)
                       X, # vector of length K
                       Y, # vector of length n (s indexed)
                       A, # vector of length K
                       Time, # scalar
                       a, # vector of length n (s indexed)
                       r, # n x I matrix (row = each s, col = each i), rate of each i in each s, 0 means that i is not in that s
                       e, # vector of length n (s indexed)
                       m # K x n matrix of {0, 1}, indicating which segment is in which background region
) { 
  K <- length(xi)
  I <- length(Lambda)
  n <- length(Y)
  
  # rate for each Ss and Bs
  Ss_rate <- (r %*% Lambda) * e # column vector of length n (s indexed)
  all_Ss_rate <- rep(1, K) %*% t(Ss_rate) # K x n matrix of S's rates. [:] for one s, identical on [..]
  all_B_rate <- (rep(1, K) %*% t(a)) * m * (xi %*% t(rep(1, n)))# K x n matrix of B's rates. [..] for one k, [:] for one s
  total_rate <- (all_Ss_rate + all_B_rate) * c(Time)
  
  # shape the matrix (Ys|xi_k)
  Y_matrix <- (rep(1, K) %*% t(Y)) * m
  
  # Xk rate
  back_rate <- xi * A * c(Time)
  
  # likelihood
  sum(dpois(x = X, # summing over k
            back_rate, 
            log = TRUE)) +
    sum(dpois(x = Y_matrix, total_rate, log = TRUE)) # summing over the entire matrix
}
