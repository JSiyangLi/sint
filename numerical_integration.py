import LaplacesDemon
import matrixStats
import cubature
import Rcpp
from scipy.stats import truncexpon, poisson, gamma
import numpy as np

# Reading relevant functions and simulation settings
# source(paste0(getwd(), "/kernel_functions.R"))
# source(paste0(getwd(), "/simulation_functions.R"))

# Reading the emcee posterior sample
emcee_results_list = LaplacesDemon.readRDS(file="/simple/emcee/results/emcee_result.RDS")
print(type(emcee_results_list$posterior_sample))

################
# Posterior kernel
################

def t_Xlik(x, t_xi, A, Time):
    xi = np.exp(t_xi)
    return poisson.logpmf(x, A * Time * xi)

def t_Lambda_ker(y, Lambda, t_mu, t_theta, t_pid, t_xi, A, Time, a, r, e):
    mu = np.exp(t_mu)
    theta = np.exp(t_theta)
    pid = LaplacesDemon.invlogit(t_pid)
    xi = np.exp(t_xi)

    log_likelihood = poisson.logpmf(y, (a * xi) + (r * e * Lambda) * Time)

    jacobian = np.log(pid) + LaplacesDemon.logdeltafunction(Lambda) + \
                np.log(1 - pid) + gamma.logpdf(Lambda, a=np.exp(2 * t_mu - t_theta), scale=np.exp(t_mu - t_theta))

    return log_likelihood + jacobian

def t_Lambda_gamma(y, Lambda, t_mu, t_theta, t_xi, A, Time, a, r, e):
    lshape = 2 * t_mu - t_theta
    lrate = t_mu - t_theta
    xi = np.exp(t_xi)

    lprior = gamma.logpdf(Lambda, a=np.exp(lshape), scale=np.exp(lrate))
    llike = poisson.logpmf(y, (a * xi) + (r * e * Lambda) * Time)

    return lprior + llike

#################
# Numerical integration
#################

def numerical(t_mu, t_theta, t_pid, t_xi, Lambda, x, y, A, Time, a, r, e):
    pid = LaplacesDemon.invlogit(t_pid)
    xi = np.exp(t_xi)

    o = np.array([t_Lambda_gamma(y, Lambda[Lambda != 0], t_mu, t_theta, t_xi, A, Time, a, r, e) for Lambda_i in Lambda[Lambda != 0]])
    io = matrixStats.logsumexp(o, axis=1) - np.log(o.shape[1])

    zeroinf = np.log(pid) + poisson.logpmf(y, lambda=a * xi * Time)
    nonzeroinf = np.log(1 - pid) + io

    return t_Xlik(x, t_xi, A, Time) + np.sum(np.vectorize(matrixStats.logsumexp)(zeroinf, nonzeroinf))

# Example usage
numerical(np.log(mu_star), np.log(theta_star), t_pid, np.log(xi_star),
          np.arange(0, 7, 0.1),
          x=30, y=np.arange(0, 3), A=A, Time=Time, a=a, r=R, e=e)
numerical(t_mu, t_theta, t_pid, t_xi,
          np.arange(0, 1, 1e-5),
          x=30, y=np.arange(0, 25), A=A, Time=Time, a=a, r=R, e=e)

##########################
# Analytical solution
##########################

def analytical(t_mu=None, t_theta=None, t_shape=2 * t_mu - t_theta,
               t_rate=t_mu - t_theta, t_pid=None, t_xi=None,
               D=None, n=None, A=None, Time=None, a=None, r=None, e=None):
    X = D[n + 1]
    Y = D[0:n]

    def summand(yi, k, t_shape, t_rate, t_xi, Time, a, r, e):
        lcombination = np.sum([LaplacesDemon.lfactorial(yi) - LaplacesDemon.lfactorial(k) - LaplacesDemon.lfactorial(yi - k)])
        lpowers = (yi - k) * (np.log(a) + t_xi + np.log(Time)) + k * (np.log(r) + np.log(e) + np.log(Time))
        llaplace = gamma.logpdf(np.exp(t_shape) + k, a=np.exp(t_shape), scale=-np.exp(t_shape - t_rate - np.log(r) - np.log(e) - np.log(Time)))
        return lcombination + lpowers + llaplace

    def prodant(yi, t_shape, t_rate, t_xi, t_pid, Time, a, r, e):
        pid = LaplacesDemon.invlogit(t_pid)
        lpois = -a * np.exp(t_xi) * Time - LaplacesDemon.lfactorial(yi)
        l0inf = yi * (np.log(a) + t_xi + np.log(Time)) + np.log(pid)
        sums = np.log(1 - pid) + np.exp(t_shape) * (t_rate) - gamma.logpdf(np.exp(t_shape)) + \
               LaplacesDemon.logsumexp(np.vectorize(summand)(yi, np.arange(0, yi + 1), np.full_like(t_shape, t_shape),
                                                            np.full_like(t_rate, t_rate), np.full_like(t_xi, t_xi),
                                                            np.full_like(Time, Time), np.full_like(a, a), np.full_like(r, r),
                                                            np.full_like(e, e)), axis=1)
        return lpois + LaplacesDemon.logplus(l0inf, sums)

    def prodant_vec(Y, t_shape, t_rate, t_xi, t_pid, Time, a, r, e):
        return np.vectorize(prodant)(Y, np.full_like(t_shape, t_shape), np.full_like(t_rate, t_rate),
                                     np.full_like(t_xi, t_xi), np.full_like(t_pid, t_pid), np.full_like(Time, Time),
                                     np.full_like(a, a), np.full_like(r, r), np.full_like(e, e))

    return poisson.logpmf(X, A * Time * np.exp(t_xi)) + np.sum(prodant_vec(Y, t_shape, t_rate, t_xi, t_pid, Time, a, r, e))

# Example usage
analytical(t_shape=t_shape, t_rate=t_rate, t_xi=t_xi, t_pid=t_pid,
           D=emcee_results_list$simulated_data, n=len(emcee_results_list$simulated_data) - 1,
           A=A, Time=Time, a=a, r=R, e=e)

################
# Example of numerical methods accuracy
################

# Debug testing: integrate()
def test_integral(Lambda, t_mu, t_theta, r, e, Time):
    mu = np.exp(t_mu)
    theta = np.exp(t_theta)
    shape = mu ** 2 / theta
    return np.exp((shape - 1) * np.log(Lambda) - Lambda * shape - r * e * Time * Lambda)

i = cubature.integrate(test_integral, np.array([0]), np.array([np.inf]), args=(t_mu, t_theta, R, e, Time))
print(np.log(i[0]))

# Numerical integration
def l_test_integral(t_mu, t_theta, t_pid, t_xi, Lambda, x, y, A, Time, a, r, e):
    pid = LaplacesDemon.invlogit(t_pid)
    xi = np.exp(t_xi)

    o = np.array([t_Lambda_gamma(y, Lambda[Lambda != 0], t_mu, t_theta, t_xi, A, Time, a, r, e) for Lambda_i in Lambda[Lambda != 0]])
    io = matrixStats.logsumexp(o, axis=1) - np.log(o.shape[1])

    zeroinf = np.log(pid) + poisson.logpmf(y, lambda=a * xi * Time)
    nonzeroinf = np.log(1 - pid) + io

    return np.vectorize(LaplacesDemon.logplus)(zeroinf, nonzeroinf)

l_test_integral(t_mu, t_theta, t_pid, t_xi, np.arange(0, 1e-0, 1e-7),
                x=emcee_results_list$simulated_data[len(emcee_results_list$simulated_data) - 1],
                y=emcee_results_list$simulated_data[0:len(emcee_results_list$simulated_data) - 1],
                A=A, Time=Time, a=a, r=R, e=e)

# Plotting
# plot(seq(0, 5e-4, by = 1e-8), exp(l_test_integral(seq(0, 5e-4, by = 1e-8), t_mu, t_theta, R, e, Time)))
# The peak is around 7e-5
# True value via Laplace:
# test_inteana <- function(k = 0, t_mu, t_theta, r, e, Time) {
#   mu <- exp(t_mu)
#   theta <- exp(t_theta)
#   shape <- mu^2 / theta
#
#   ifelse((shape + k) > 0, lgamma(shape) + (-shape) * log(shape + r * e * Time), -Inf)
# }
# test_inteana(t_mu = t_mu, t_theta = t_theta, r = R, e = e, Time = Time)
# LESSON: Focus on the peak and don't include delta functions (that changes the weighting of the whole integral).

##################
# Complete integrated kernel (partial marginal likelihood)
##################
# Remaining code for the integrated kernel is not provided.
