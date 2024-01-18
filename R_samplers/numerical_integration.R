library(LaplacesDemon) # for logistic function
library(matrixStats)
library(cubature)
library(Rcpp)
sourceCpp("logsum.cpp")

# reading relevant functions and simulation settings
#source(paste0(getwd(), "/kernel_functions.R"))
#source(paste0(getwd(), "/simulation_functions.R"))

# reading the emcee posterior sample
emcee_results_list <- readRDS(file = paste0(getwd(), "/simple/emcee/results/emcee_result.RDS"))
class(emcee_results_list$posterior_sample)

################
# posterior kernel
###############
t_Xlik <- function(x, t_xi, A, Time) {
  xi <- exp(t_xi)
  dpois(x = x,
        A * Time * xi, 
        log = TRUE)
}

t_Lambda_ker <- function(y, Lambda, t_mu, t_theta, t_pid, t_xi,
                         A, Time, a, r, e) {
  #Lambda <- exp(t_Lambda)
  mu <- exp(t_mu)
  theta <- exp(t_theta)
  pid <- invlogit(t_pid)
  xi <- exp(t_xi)
  
  dpois(x = y, 
        (c(a * xi) + c(r) * c(e) * Lambda) * c(Time), 
        log = TRUE) + #Jacobian:
    #t_Lambda + # because Lambda is fed into the function, not log Lambda, so the Jacobian is not needed.
    sapply(Lambda, 
           function(i) logplus(log(pid) + logdeltafunction(i),
                               log(1 - pid) + dgamma(x = i,
                                                     shape = exp(2 * t_mu - t_theta),
                                                     rate = exp(t_mu - t_theta),
                                                     log = TRUE)))
}
# plotting
#plot(seq(0, 5e-3, by = 1e-7), exp(t_Lambda_ker(200, seq(0, 5e-3, by = 1e-7), t_mu, t_theta, t_pid, t_xi, A, Time, a, R, e))) # peak depends on data.

t_Lambda_gamma <- function(y, Lambda, t_mu, t_theta, t_xi, A, Time, a, r, e) {
  lshape <- 2 * t_mu - t_theta
  lrate <- t_mu - t_theta
  xi <- exp(t_xi)
  
  lprior <- dgamma(Lambda, shape = exp(lshape), rate = exp(lrate), log = TRUE)
  llike <- dpois(x = y, 
                 (c(a * xi) + c(r) * c(e) * Lambda) * c(Time), 
                 log = TRUE)
  lprior + llike
}

#########
# numerical integration
#########
numerical <- function(t_mu, t_theta, t_pid, t_xi, Lambda, # Lambda > 0 needed
                      x, y, A, Time, a, r, e) {
  pid <- invlogit(t_pid)
  xi <- exp(t_xi)
  
  o <- outer(y, Lambda[Lambda != 0], FUN = t_Lambda_gamma,
             t_mu = t_mu, t_theta = t_theta, t_xi = t_xi, A = A, Time = Time, a = a, r = r, e = e)
  io <- apply(o, 1, logplusvec) - log(ncol(o))
  zeroinf <- log(pid) + dpois(y, lambda = a * xi * Time, log = TRUE)
  nonzeroinf <- log(1 - pid) + io
  t_Xlik(x, t_xi, A, Time) + sum(sapply(1:length(zeroinf), function(i) logplus(zeroinf[i], nonzeroinf[i])))
  
  
  # test
  #list(t_Xlik(x, t_xi, A, Time),
  #     io)
}


##########################
# analytical solution
#########################
analytical <- function(t_mu = NULL, t_theta = NULL, 
                       t_shape = 2 * t_mu - t_theta, 
                       t_rate = t_mu - t_theta, t_pid, t_xi, 
                       D, n = length(D) - 1, A, Time, a, r, e) {
  # reshape formal arguments
  X <- D[n + 1]
  Y <- D[1:n]
  
  # construction
  summand <- function(yi, k, t_shape, t_rate, t_xi, Time,
                      a, r, e) {
    lcombination <- lfactorial(yi) - (lfactorial(k) + lfactorial(yi-k))
    lpowers <- (yi - k) * (log(a) + t_xi + log(Time)) + k * (log(r) + log(e) + log(Time))
    llaplace <- lgamma(exp(t_shape) + k) + (-exp(t_shape) - k) * (logplus(t_rate, log(r) + log(e) + log(Time)))
    lcombination + lpowers + ifelse((exp(t_shape) + k) > 0, llaplace, -Inf) # the condition of Laplace transform
  }
  prodant <- function(yi, t_shape, t_rate, t_xi, t_pid, Time,
                      a, r, e) {
    pid <- invlogit(t_pid)
    
    lpois <- -a * exp(t_xi) * Time - lfactorial(yi)
    l0inf <- yi * (log(a) + t_xi + log(Time)) + log(pid)
    sums <- log(1 - pid) + exp(t_shape) * (t_rate) - lgamma(exp(t_shape)) + 
      logplusvec(summand(yi = yi, k = 0:yi, t_shape = t_shape, t_rate = t_rate, t_xi = t_xi, Time = Time, a = a, r = r, e = e))
    lpois + logplus(l0inf, sums)
    
    ############## test
    #dpois(yi, lambda = a * exp(t_xi) * Time, log = TRUE) + log(pid)
  }
  prodant_vec <- function(Y, t_shape, t_rate, t_xi, t_pid, Time, 
                          a, r, e) {
    sapply(Y, prodant, t_shape = c(t_shape), t_rate = c(t_rate), t_pid = c(t_pid), t_xi = c(t_xi), Time = c(Time), a = c(a), r = c(r), e = c(e))
  }
  dpois(x = X, lambda = A * Time * exp(t_xi), log = TRUE) + 
    sum(prodant_vec(Y = Y, t_shape = t_shape, t_rate = t_rate, t_pid = t_pid, t_xi = t_xi, Time = Time, a = a, r = r, e = e))
  # test
  #list(dpois(x = X, lambda = A * Time * exp(t_xi), log = TRUE),
  #     prodant_vec(Y = Y, t_mu = t_mu, t_theta = t_theta, t_pid = t_pid, t_xi = t_xi, Time = Time, a = a, r = r, e = e))
}


##################
# complete integrated kernel (partial marginal likelihood)
#################