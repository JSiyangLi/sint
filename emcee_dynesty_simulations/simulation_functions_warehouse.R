#######################
# broken power law lambdas
rbrokenpower <- function(n, pllocation_vec, plshape_vec = rep(1, m),
                         m = length(pllocation_vec)) {
  # a single broken power law distribution random number
  rmixtruncpareto_single <- function(pllocation_vec, plshape_vec = rep(1, m),
                                     m = length(pllocation_vec)) {
    if (length(plshape_vec) != m)
      stop("The length of shape parameters do not match the number of components")
    
    pltau_diff <- -plshape_vec[-m] * diff(log(pllocation_vec)) # Irina thesis equation 2.7
    initial_logpjs <- sapply(pltau_diff, logminus, x = 0) + cumsum(pltau_diff)
    last_logpj <- logminus(0, logplusvec(initial_logpjs))
    mix_weights <- exp(c(initial_logpjs, last_logpj)) / sum(exp(c(initial_logpjs, last_logpj))) # normalize to even out computation error
    mix_truncparetos <- rtruncpareto(n = m - 1, lower = pllocation_vec[-m], upper = pllocation_vec[-1], shape = plshape_vec[-m])
    mix_pareto <- VGAM::rpareto(n = 1, scale = pllocation_vec[m], shape = plshape_vec[m])
    mix_weights %*% c(mix_truncparetos, mix_pareto)
  }
  
  replicate(n = n, expr = rmixtruncpareto_single(pllocation_vec, plshape_vec, m))
}

Lambdastar_simulator_brokenpower <- function(pid, pllocation_vec, plshape_vec) {
  ########
  # step 2
  ########
  zeroindexes <- which(log(runif(n)) < log(pid))
  Lambda_stars <- rbrokenpower(n, pllocation_vec, plshape_vec)
  Lambda_stars[zeroindexes] <- 0
  
  Lambda_stars
}

###################
# prior sampling functions
###################
source("kernel_functions.R")
prior_sampling <- function(sample_size = 1, 
                           muloc = 1.2e-6, thetaloc = 8.7e-12, 
                           shape1 = 1, shape2 = 1,
                           mu0 = 10^6 / (A * Time), theta0 = 10^18 / (A * Time)^2) {
  # truncated Cauchy from Programs in R for Computing Truncated Cauchy Distributions, Saralees Nadarajah1 and Samuel Kotz (2007)
  qtcauchy <- function(p, location = 0, scale = 1, a, b) {
    second <- ifelse(missing(a), 0, pcauchy(q = a,location = location,scale = scale))
    first <- ifelse(missing(b), 1, pcauchy(q = b,location = location,scale = scale))
    qcauchy(second + p * (first - second),location = location,scale = scale)
  }
  rtcauchy = function(n,location=0,scale=1,a,b){qtcauchy(p = runif(n,min=0,max=1),location,scale,a,b)}
  
  # mu
  mu_sample <- rtcauchy(n = sample_size, location = muloc, scale = 100 * muloc, a = 0)
  # theta
  theta_sample <- rtcauchy(n = sample_size, location = thetaloc, scale = 100 * thetaloc, a = 0)
  # pid
  pid_sample <- rbeta(n = sample_size, shape1 = shape1, shape2 = shape2)
  # xi # too hard to sample from tiny gamma
  uniflow <- 1
  while (uniflow >= 1e-300) {
    tgamma <- rgamma(n = min(1e6, sample_size * 1e3), shape = (mu0^2) / theta0, rate = mu0 / theta0)
    uniflow <- max(min(tgamma[tgamma != 0]), .Machine$double.xmin) # giving the minimum probability value in the sampling.
  } # this restricts the smallest rgamma value from being too small, hence stops them from being 0.
  plowbound <- pgamma(q = uniflow, shape = (mu0^2) / theta0, rate = mu0 / theta0, log.p = TRUE)
  xi_sample <- qgamma(p = log(runif(n = sample_size, min = exp(plowbound), max = 1)), 
                      shape = (mu0^2) / theta0, rate = mu0 / theta0, log.p = TRUE)
  
  cbind(mu_sample, theta_sample, pid_sample, xi_sample)
}
# test
prior_sampling(sample_size = 10)

# broken power law
################################
# tests from numerical_integration
###########
numerical(log(mu_star), log(theta_star), t_pid, log(xi_star), #t_Lambda = -Inf,  
          #t_Lambda = seq(-31, 31, by = 0.01), 
          Lambda = seq(0, exp(7), by = 0.1),
          #Lambda = 0,
          x = 30, y = 0:3, A, Time, a, R, e = e)
numerical(t_mu, t_theta, t_pid, t_xi, #t_Lambda = -Inf,  
          #t_Lambda = seq(-31, 31, by = 0.01), 
          Lambda = seq(0, 1, by = 1e-5),
          #Lambda = 0,
          x = 30, y = 0:25, A, Time, a, R, e = e)

# example
analytical(D = c(0:36, 250000), t_shape = t_shape, t_rate = t_rate, t_xi = t_xi, t_pid = t_pid, A = A, Time = Time, a = a, r = R, e = e)
################
# example of numerical methods accuracy
#################
# debug testing: integrate():
test_integral <- function(Lambda, t_mu, t_theta, r, e, Time) {
  mu <- exp(t_mu)
  theta <- exp(t_theta)
  shape <- mu^2 / theta
  
  exp((shape - 1) * log(Lambda) - Lambda * shape - r * e * Time * Lambda)
}
(i <- integrate(test_integral, lower = 0, upper = Inf, t_mu = t_mu, t_theta = t_theta, r = R, e = e, Time = Time))
log(i$value)

# numerical integration:
l_test_integral <- function(t_mu, t_theta, t_pid, t_xi, Lambda, # Lambda > 0 needed
                            x, y, A, Time, a, r, e) {
  pid <- invLaplacesDemon::logit(t_pid)
  xi <- exp(t_xi)
  
  o <- outer(y, Lambda[Lambda != 0], FUN = t_Lambda_gamma,
             t_mu = t_mu, t_theta = t_theta, t_xi = t_xi, A = A, Time = Time, a = a, r = R, e = e)
  io <- apply(o, 1, logplusvec) - log(ncol(o))
  zeroinf <- log(pid) + dpois(y, lambda = a * xi * Time, log = TRUE)
  nonzeroinf <- log(1 - pid) + io
  sapply(1:length(zeroinf), function(i) logplus(zeroinf[i], nonzeroinf[i]))
  
  
  # test
  #list(t_Xlik(x, t_xi, A, Time),
  #     io)
}
#l_test_integral(t_mu, t_theta, t_pid, t_xi, seq(0, 1e-0, by = 1e-7), 
#                x = emcee_results_list$simulated_data[length(emcee_results_list$simulated_data)], y = emcee_results_list$simulated_data[-length(emcee_results_list$simulated_data)], 
#                A = A, Time = Time, a = a, r = R, e = e)
#plot(seq(0, 5e-4, by = 1e-8), exp(l_test_integral(seq(0, 5e-4, by = 1e-8), t_mu, t_theta, R, e, Time))) # the peak is around 7e-5
# true value via Laplace:
#test_inteana <- function(k = 0, t_mu, t_theta, r, e, Time) {
#  mu <- exp(t_mu)
#  theta <- exp(t_theta)
#  shape <- mu^2 / theta

#  ifelse((shape + k) > 0, lgamma(shape) + (-shape) * log(shape + r * e * Time), -Inf)
#}
#test_inteana(t_mu = t_mu, t_theta = t_theta, r = R, e = e, Time = Time)
# LESSON: focus on the peak and don't include delta functions (that changes the weightening of the whole integral).


##############
# tests from kernel_functions
################
integrated_kernel(t_shape = t_shape, t_rate = t_rate, t_pid = t_pid, t_xi = t_xi, D = c(0:3, 30), A = A, Time = Time, a = a, r = R, e = e)

