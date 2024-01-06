import numpy as np
from scipy.stats import poisson, dgamma, cauchy, beta
from scipy.special import gammaln, expit, log1p
import ctypes

####
# load related files
####
# in terminal, run cd /Users/jl1023/Library/CloudStorage/OneDrive-ImperialCollegeLondon/PhD/source_intensities_git
# then run g++ -fPIC -shared -o logsum.so logsum.cpp
from simulation_functions import t_mu, t_theta, A, Time, a, R, e
from numerical_integration import analytical, numerical
clibrary = ctypes.CDLL("/logsum.so")
def logplus(x, y):
    return clibrary.py_logplus(x, y)
def logplusvec(x):
    return clibrary.py_logplusvec(x)
def logminus(x, y):
    return clibrary.py_logminus(x, y)


#####################
# Truncated Pareto density function
#####################

def dtruncpareto(x, lower, upper, shape, log=False):
    """
    Truncated Pareto density function.

    Parameters:
    - x: array-like, values at which to evaluate the density function
    - lower: array-like, lower bound of the truncated Pareto distribution
    - upper: array-like, upper bound of the truncated Pareto distribution
    - shape: array-like, shape parameter of the truncated Pareto distribution
    - log: bool, if True, return the log-density; otherwise, return the density

    Returns:
    - array-like, values of the density function at the given points
    """
    if not isinstance(log, bool):
        raise ValueError("bad input for argument 'log'")

    L = max(len(x), len(lower), len(upper), len(shape))
    x = np.broadcast_to(x, L)
    lower = np.broadcast_to(lower, L)
    upper = np.broadcast_to(upper, L)
    shape = np.broadcast_to(shape, L)

    logdensity = np.full_like(x, np.nan, dtype=float)
    xok = (0 < lower) & (lower <= x) & (x <= upper) & (shape > 0)

    logdensity[xok] = np.log(shape[xok]) + \
                      shape[xok] * np.log(lower[xok]) - \
                      (shape[xok] + 1) * np.log(x[xok]) - \
                      np.log1p(-np.power(lower[xok] / upper[xok], shape[xok]))

    logdensity[shape <= 0] = np.nan
    logdensity[upper < lower] = np.nan
    logdensity[0 > lower] = np.nan

    return logdensity if log else np.exp(logdensity)


#################################################################
# Posterior kernel
#################################################################

#############
# Likelihood
#############

def loglik(xi, Lambda, D, A, Time, a, r, e):
    """
    Log-likelihood function.

    Parameters:
    - xi: float, parameter
    - Lambda: float, parameter
    - D: array-like, data
    - A: float, parameter
    - Time: float, parameter
    - a, r, e: array-like, parameters

    Returns:
    - float, log-likelihood value
    """
    if xi >= np.finfo(float).tiny:
        return poisson.logpmf(D[-1], A * Time * xi) + \
               np.sum(poisson.logpmf(D[:-1], (a * xi) + (r * e * Lambda) * Time))
    else:
        return -1e100


def t_loglik(t_xi, t_Lambda, D, A, Time, a, r, e):
    xi = np.exp(t_xi)
    Lambda = np.exp(t_Lambda)
    return poisson.logpmf(D[-1], A * Time * xi) + \
           np.sum(poisson.logpmf(D[:-1], (a * xi) + (r * e * Lambda) * Time))


def loglik_ext(xi, Lambda, X, Y, A, Time, a, r, e, m):
    """
    Log-likelihood function for multiple combinations.

    Parameters:
    - xi: array-like, vector of length K
    - Lambda: array-like, vector of length n (i indexed)
    - X: array-like, vector of length K
    - Y: array-like, vector of length n (s indexed)
    - A: array-like, vector of length K
    - Time: float, scalar
    - a: array-like, vector of length n (s indexed)
    - r: array-like, n x n matrix (row = each s, col = each i), rate of each i in each s, 0 means that i is not involved in that s
    - e: array-like, vector of length n (s indexed)
    - m: array-like, K x n matrix of {0, 1}, indicating which segment is in which region

    Returns:
    - float, log-likelihood value
    """
    return np.sum(poisson.logpmf(X, A * Time * xi)) + \
           np.sum(poisson.logpmf(Y, (a * xi) + (r * e * Lambda) * Time))

#####################
# Complete integrated kernel
########################

def integrated_kernel(D, A, Time, a, r, e, t_mu=None, t_theta=None, t_shape=2 * t_mu - t_theta,
                      t_rate=t_mu - t_theta, t_pid=None, t_xi=None,
                      shape1=1, shape2=1, muloc=1.2e-6, thetaloc=8.7e-12,
                      mu0=10**6 / (A * Time), theta0=10**18 / (A * Time)**2,
                      analytical_method=True, Lambda_seq=None):
    """
    Complete integrated kernel function.

    Parameters:
    - t_mu, t_theta, t_shape, t_rate, t_pid, t_xi: parameters
    - D: array-like, data
    - n: int, length of data minus 1
    - A, Time: float, parameters
    - a, r, e: array-like, parameters
    - shape1, shape2: parameters
    - muloc, thetaloc: parameters
    - mu0, theta0: parameters
    - analytical_method: bool, use analytical or numerical method
    - Lambda_seq: array-like, sequence

    Returns:
    - float, result of the integrated kernel
    """
    n = len(D) - 1

    integral = analytical(t_shape=t_shape, t_rate=t_rate, t_xi=t_xi, t_pid=t_pid, D=D,
                          A=A, Time=Time, a=a, r=r, e=e) if analytical_method else \
               numerical(t_mu = t_mu, t_theta = t_theta, t_pid = t_pid, t_xi = t_xi, Lambda=Lambda_seq, x=D[n + 1], y=D[1:n + 1],
                         A=A, Time=Time, a=a, r=r, e=e)

    hyperpriors = t_logpidhypprior(t_pid, shape1, shape2) + \
                  t_logshaperatehypprior(t_shape, t_rate, muloc, thetaloc) + \
                  t_logxiprior_gamma(t_xi, mu0, theta0)

    if t_mu is not None and t_theta is not None:
        hyperpriors = t_logpidhypprior(t_pid, shape1, shape2) + \
                      t_logmuhypprior(t_mu, muloc) + \
                      t_logthetahypprior(t_theta, thetaloc) + \
                      t_logxiprior_gamma(t_xi, mu0, theta0)

    return integral + hyperpriors

import numpy as np

def t_int_kernel_gamma(obj, D, A, Time, a, R, e, mu_theta_param=False):
    # vectorising a function to make it 'S3'
    return np.vectorize(lambda zeta: t_int_kernel_gamma.numeric(zeta, D, A, Time, a, R, e, mu_theta_param))
def t_int_kernel_gamma_numeric(zeta, D, A, Time, a, R, e, mu_theta_param=False):
    if mu_theta_param:
        ll = integrated_kernel(t_mu=zeta[0],
                               t_theta=zeta[1],
                               t_pid=zeta[2],
                               t_xi=zeta[3],
                               D=D,
                               A=A,
                               Time=Time,
                               a=a,
                               r=R,
                               e=e)
    else:
        ll = integrated_kernel(t_shape=zeta[0],
                               t_rate=zeta[1],
                               t_pid=zeta[2],
                               t_xi=zeta[3],
                               D=D,
                               A=A,
                               Time=Time,
                               a=a,
                               r=R,
                               e=e)

    return -np.inf if np.isnan(ll) or np.isinf(ll) else ll

############
# prior
# gamma
#############

def logdeltafunction(i):
    return 0 if i == 0 else -np.inf

def logLambdaprior_gamma(Lambda, mu, theta, pid):
    return np.sum([np.log(logplus(np.log(pid) + logdeltafunction(i),
                                  np.log(1 - pid) + dgamma(x=i, shape=mu**2 / theta, rate=mu / theta, log=True)))
                   for i in Lambda])

def logxiprior_gamma(xi, mu0, theta0):
    return dgamma(xi, shape=mu0**2 / theta0, rate=mu0 / theta0, log=True)



####################
# Prior transform
####################

def prior_transform(u):
    """
    Prior transform function.

    Parameters:
    - u: array-like, parameters

    Returns:
    - array-like, transformed parameters
    """
    if np.any(u < np.finfo(float).tiny) or np.any(u >= 1):
        return -np.inf
    else:
        muloc = 1.2e-6
        thetaloc = 8.7e-12
        lmu0 = np.log(10**6) - np.log(A) - np.log(Time)
        ltheta0 = np.log(10**18) - 2 * (np.log(A) + np.log(Time))

        umu = truncexpon.ppf(u[0], loc=muloc, scale=100 * muloc)
        utheta = truncexpon.ppf(u[1], loc=thetaloc, scale=100 * thetaloc)
        upid = u[2]
        uxi = qgamma_bounded(u[3], shape=np.exp(2 * lmu0 - ltheta0), rate=np.exp(lmu0 - ltheta0))

        return np.array([umu, utheta, upid, uxi, 0, 0])

####################
# Posterior kernel
####################
def logpostker_unified_gamma(params, D, A, Time, a, r, e):
    """
    Unified posterior kernel function.

    Parameters:
    - params: array-like, parameters
    - D: array-like, data
    - A, Time: float, parameters
    - a, r, e: array-like, parameters

    Returns:
    - float, result of the unified posterior kernel
    """
    return integrated_kernel(t_mu=params[0], t_theta=params[1], t_shape=params[2], t_rate=params[3],
                             t_pid=params[4], t_xi=params[5], D=D, A=A, Time=Time, a=a, r=R, e=e,
                             shape1=1, shape2=1, muloc=1.2e-6, thetaloc=8.7e-12,
                             mu0=10**6 / (A * Time), theta0=10**18 / (A * Time)**2)


def t_logpostker_unified_gamma(params, D, A, Time, a, r, e):
    """
    Transformed unified posterior kernel function.

    Parameters:
    - params: array-like, parameters
    - D: array-like, data
    - A, Time: float, parameters
    - a, r, e: array-like, parameters

    Returns:
    - float, result of the transformed unified posterior kernel
    """
    return logpostker_unified_gamma(*params, D=D, A=A, Time=Time, a=a, r=R, e=e)
