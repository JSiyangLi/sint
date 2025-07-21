########################
### FUNCTIONS        ###
########################

# NEGATIVE Conditional Dist'n of LOG mean, LOG var, LOGIT pi of zero-inflated 
# gamma pop dist'n
# This version is used in the Partially Collapsed Gibbs .
# It uses likelihood based on zero-inflated neg binomial dist'n of (missing)
# source counts, after marginalizing out lambda = Poisson intensity parameters
# and Cauchy priors on gamma mean and variance given in Lazhi's paper.
#
pcg.cond.post <- function(pars, cnts, r){
  # pars = (log(mu), log(theta), logit(pi)) in Lazhi's notation and priors
  # cnts = zero-inflated negative binomial data (counts)
  # r    = exposures, i.e., conditional on lambda, data ~ Pois(lambda *r)
  
  # LOG Cauchy priors on mu and theta
  prior <- dcauchy(exp(pars[1]), location=1.2e-06, scale=1.2e-04, log=TRUE)
  prior <- prior + dcauchy(exp(pars[2]), location=8.7e-12, scale=8.7e-10, log=TRUE)
  prior <- prior + pars[1] +pars[2]   # LOG Jacobian of log transforms
  
  # Uniform prior on pi, LOG Jacobian on transform to logit(pi)
  prior <- prior + pars[3] - 2*log(1 + exp(pars[3])) 
  
  # Likelihood
  alpha <- exp(2*pars[1]-pars[2])
  beta  <- exp(pars[1]-pars[2])
  pi    <- exp(pars[3]) / (1+ exp(pars[3]))
  loglike <-  log(1-pi) + lgamma(alpha+cnts) -lgamma(alpha) - lgamma(cnts+1) +
      alpha * log(beta) + cnts*log(r) - (alpha+cnts)*log(r+beta)
  loglike[cnts==0] <- log(exp(loglike[cnts==0]) + pi) 
    
  -(sum(loglike) +prior)
}


# NEGATIVE Conditional Dist'n of mean and var (mu, theta) of gamma pop dist'n
# This version is used in the Full Gibbs Sampler.
# It uses likelihood based on gamma dist'n of Poisson intensity parameters
# and Cauchy priors on gamma mean and variance given in Lazhi's paper.
#
cond.post <- function(pars, data){
  #pars = (mu, theta) in Lazhi's notation and priors
  value <- -Inf
  if((pars[1] >0)*(pars[2]>0)){
    value <- dcauchy(pars[1], location=1.2e-06, scale=1.2e-04, log=TRUE)
    value <- value + dcauchy(pars[2], location=8.7e-12, scale=8.7e-10, log=TRUE)
    value <- value +
      sum(dgamma(data, pars[1]^2/pars[2], pars[1]/pars[2], log=TRUE))
  }
  -value
}

# NEGATIVE Conditional Dist'n of LOG mean and LOG var of gamma pop dist'n
# This version is used in the Full Gibbs Sampler.
# It uses likelihood based on gamma dist'n of Poisson intensity parameters
# and Cauchy priors on gamma mean and variance given in Lazhi's paper.
# The only difference between this and cond.post is that mu and theta are 
# transfomred to log(mu) and log(theta)
#
cond.post.log <- function(pars, data){
  #pars = (log(mu), log(theta))  in Lazhi's notation and priors
  value <- dcauchy(exp(pars[1]), location=1.2e-06, scale=1.2e-04, log=TRUE)
  value <- value + dcauchy(exp(pars[2]), location=8.7e-12, scale=8.7e-10, log=TRUE)
  value <- value + pars[1] +pars[2]   # Jacobian 
  value <- value +
    sum(dgamma(data, exp(2*pars[1]-pars[2]), exp(pars[1]-pars[2]), log=TRUE))
  return(-value)
}

# A version of Gibbs step that uses the log population parameters
#for(iter in 1:1000){
  # convert current shape and rate to log(mean) and log(variance)
#  current <- c(log(params$pop.shape[iter]) - log(params$pop.rate[iter]), 
#               log(params$pop.shape[iter]) - 2*log(params$pop.rate[iter]))
#  fit <- optim(current, cond.post, data=not.zero.rates, hessian=TRUE)
#  variance <- solve(fit$hessian)
#  proposal <- rmvnorm( 1, mean=fit$par, sigma=variance )
#  accept.prob <- exp( -cond.post(proposal, data=not.zero.rates) 
#                      + dmvnorm(current, mean=fit$par, sigma=variance, log=TRUE ) 
#                      + cond.post(current, data=not.zero.rates)  
#                      - dmvnorm(proposal, mean=fit$par, sigma=variance, log=TRUE ))
#  print(c(current, proposal, accept.prob))
#  params$pop.shape[iter+1] <- params$pop.shape[iter]
#  params$pop.rate[iter+1]  <- params$pop.rate[iter]
#  if(runif(1) < accept.prob){
#    params$pop.shape[iter+1] <- exp(2*proposal[1] - proposal[2]) 
#    params$pop.rate[iter+1] <- exp(proposal[1] - proposal[2])
#  }
#}