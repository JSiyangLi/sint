from scipy.stats import poisson, gamma, cauchy, beta, expon
from scipy.special import gammaln, expit, log1p, gammaincinv, logit
import math
import sys
import numpy as np
import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
utils = importr('utils')
base = importr('base')
stats = importr('stats')
LaplacesDemon = importr('LaplacesDemon')
#from rpy2.robjects import numpy2ri
#from rpy2.robjects import default_converter
#np_cv_rules = default_converter + numpy2ri.converter

####
# load related files
####
from simulation_functions import t_mu, t_theta, A, Time, a, R, e
from numerical_integration import analytical, numerical, logplus, logplusvec, logminus, logdeltafunction
import simulation_functions
#############
# Likelihood
#############

def loglik(xi, Lambda, D, A = None, Time = None, a = None, r = None, e = None):  # pass
    if A is None:
        A = simulation_functions.A
    if Time is None:
        Time = simulation_functions.Time
    if a is None:
        a = simulation_functions.a
    if r is None:
        r = simulation_functions.R
    if e is None:
        e = simulation_functions.e
    if xi >= np.finfo(float).eps:
        if isinstance(Lambda, (np.floating, np.integer, float, int)):
            log_likelihood = poisson.logpmf(D[-1], A * Time * xi) + np.sum(poisson.logpmf(D[:-1], (a * xi + r * e * Lambda) * Time))
        elif isinstance(Lambda, (list, np.ndarray)):
            log_likelihood = poisson.logpmf(D[-1], A * Time * xi) + np.sum(poisson.logpmf(D[:-1], [(a * xi + r * e * Lambda_j) * Time for Lambda_j in Lambda]))
        else:
            log_likelihood = np.nan
        return log_likelihood
    else:
        return -1e100


def t_loglik(t_xi, t_Lambda, D, A, Time, a, r, e): # pass
    xi = np.exp(t_xi)

    Lambda = np.exp(t_Lambda)
    if isinstance(Lambda, (np.floating, np.integer, float, int)):
        log_likelihood = poisson.logpmf(D[-1], A * Time * xi) + np.sum(
            poisson.logpmf(D[:-1], (a * xi + r * e * Lambda) * Time))
    elif isinstance(Lambda, (list, np.ndarray)):
        log_likelihood = poisson.logpmf(D[-1], A * Time * xi) + np.sum(
            poisson.logpmf(D[:-1], [(a * xi + r * e * Lambda_j) * Time for Lambda_j in Lambda]))
    else:
        log_likelihood = np.nan
    return log_likelihood

def loglik_ext(xi, Lambda, X, Y, A, Time, a, r, e, m): # not needed yet
    K, n = len(xi), len(Lambda)
    region_likelihoods = np.zeros((n, K))

    for i in range(K):
        region_likelihoods[:, i] = poisson.logpmf(X, A[i] * Time * np.exp(xi[i]))

    for s in range(n):
        region_likelihoods[s, :] += np.sum(poisson.logpmf(Y, (a[s] * np.exp(xi) + np.dot(r[s, :], e * Lambda)) * Time) * m[:, s])

    return np.sum(region_likelihoods)

#####################
# Complete integrated kernel
########################

def integrated_kernel(D, A, Time, a, r, e, t_mu, t_theta, t_pid, t_xi,
                      analytical_method=True, Lambda_seq=None): # pass
    shape1 = 1
    shape2 = 1
    mu0 = 10 ** 6 / (A * Time)
    theta0 = 10 ** 18 / (A * Time) ** 2
    n = len(D) - 1

    integral = analytical(t_mu = t_mu, t_theta = t_theta, t_xi=t_xi, t_pid=t_pid, D=D,
                          A=A, Time=Time, a=a, r=r, e=e) if analytical_method else \
        numerical(t_mu=t_mu, t_theta=t_theta, t_pid=t_pid, t_xi=t_xi, Lambda=Lambda_seq, x=D[n + 1], y=D[1:n + 1],
                  A=A, Time=Time, a=a, r=r, e=e)

    # Other functions
    hyperpriors = t_logpidhypprior(t_pid, shape1, shape2) + t_logmuhypprior(t_mu) + t_logthetahypprior(t_theta) + t_logxiprior_gamma(t_xi, mu0, theta0)

    return integral + hyperpriors

def t_int_kernel_gamma(zeta, D, A = None, Time = None, a = None, R = None, e = None): # pass
    if A is None:
        A = simulation_functions.A
    if Time is None:
        Time = simulation_functions.Time
    if a is None:
        a = simulation_functions.a
    if R is None:
        R = simulation_functions.R
    if e is None:
        e = simulation_functions.e
    def t_int_kernel(zeta):
        return integrated_kernel(t_mu=zeta[0],
                                 t_theta=zeta[1],
                                 t_pid=zeta[2],
                                 t_xi=zeta[3],
                                 D=D,
                                 A=A,
                                 Time=Time,
                                 a=a,
                                 r=R,
                                 e=e)

    if np.ndim(zeta) == 1:
        return [t_int_kernel(zeta)]
    elif np.ndim(zeta) == 2:
        return np.apply_along_axis(t_int_kernel, axis=1, arr=zeta)

############
# prior
# gamma
#############

def logLambdaprior_gamma(Lambda, mu, theta, pid): # pass
    return np.sum([logplus(np.log(pid) + logdeltafunction(i),
                                  np.log(1 - pid) + gamma.logpdf(x=i, a=mu**2 / theta, scale=theta / mu))
                   for i in Lambda])

def logxiprior_gamma(xi, mu0, theta0): # pass
    return gamma.logpdf(xi, a=mu0**2 / theta0, scale=theta0 / mu0)

def t_logLambdaprior_gamma(t_Lambda, t_mu, t_theta, t_pid): # pass
    Lambda = np.exp(t_Lambda)
    pid = 1 / (1 + np.exp(-t_pid))

    def zeroinfgamma(i):
        output = logplus(logdeltafunction(i) + np.log(pid),
                         np.log(1 - pid) + gamma.logpdf(x=i,
                                                        a=np.exp(2 * t_mu - t_theta),
                                                        scale=np.exp(t_theta - t_mu)))
        return output

    log_prior_slot = [zeroinfgamma(j) for j in Lambda]
    log_prior = np.sum(t_Lambda)
    log_prior += np.sum(log_prior_slot)

    # Jacobian: prod(Lambda) => log Jacobian: sum(t_Lambda)
    #log_prior += np.sum(t_Lambda)

    return log_prior

def t_logxiprior_gamma(t_xi, mu0, theta0): # pass
    xi = np.exp(t_xi)
    log_prior = np.array(gamma.logpdf(x = xi,
                                      a=mu0**2 / theta0,
                                      scale=theta0 / mu0)) + t_xi

    return log_prior

############
# hyperprior
############
def logpidhypprior(pid, shape1=1, shape2=1): # pass
    return beta.logpdf(pid, a=shape1, b=shape2)

def logmuhypprior(mu, muloc): # pass
    if mu < 0:
        return float('-inf')
    return cauchy.logpdf(mu, loc=muloc, scale=100 * muloc)

def logthetahypprior(theta, thetaloc): # pass
    if theta < 0:
        return float('-inf')
    return cauchy.logpdf(theta, loc=thetaloc, scale=100 * thetaloc)

def t_logpidhypprior(t_pid, shape1=1, shape2=1): # pass
    pid = 1 / (1 + np.exp(-t_pid))
    return beta.logpdf(pid, a=shape1, b=shape2) + np.log(pid) + np.logaddexp(0, np.log(pid) - t_pid)

def t_logmuhypprior(t_mu, muloc=1.2e-6): # pass
    mu = np.exp(t_mu)
    return cauchy.logpdf(mu, loc=muloc, scale=100 * muloc) + t_mu

def t_logthetahypprior(t_theta, thetaloc=8.7e-12): # pass
    theta = np.exp(t_theta)
    return cauchy.logpdf(theta, loc=thetaloc, scale=100 * thetaloc) + t_theta

def t_logshaperatehypprior(t_shape, t_rate, muloc=1.2e-6, thetaloc=8.7e-12): # pass
    shapeloc = muloc**2 / thetaloc
    scaleloc = muloc / thetaloc

    return (cauchy.logpdf(np.exp(t_shape), loc=shapeloc, scale=100 * shapeloc) +
            cauchy.logpdf(np.exp(t_rate), loc=scaleloc, scale=100 * scaleloc) +
            t_shape + t_rate)

#############
# joint prior
#############
def q0inf_gamma(prob_vec, shape, rate, pid): # pass
    output_vec = np.full_like(prob_vec, np.nan)
    output_vec[prob_vec <= pid] = 0

    gamma_part_p = prob_vec[prob_vec > pid] - pid
    trans_p = gammaln(shape) + np.log(gamma_part_p) - np.log(1 - pid)
    #print(shape)
    #print("gammaln", gammaln(shape), "gammapart", np.log(gamma_part_p), "nplog", np.log(1 - pid))
    output = gammaincinv(shape, np.exp(trans_p)) / rate
    output_vec[prob_vec > pid] = output

    return output_vec

# Function to compute qgamma_bounded
def qgamma_bounded(p, shape, rate): # pass
    plowbound = gamma.cdf(np.finfo(float).eps, a=shape, scale=1/rate)
    # transform the input p so it is bounded below by machine xmin
    pt = (1 - plowbound) * p + plowbound
    return gamma.ppf(pt, a=shape, scale=1/rate)

# Function to sample from gamma distributions with tiny shape and rate on log scale
def rgamma_transformed(n, shape, rate): # pass
    # Uniform to Exponential: gexp ~ exp(shape)
    gexp = expon.ppf(np.random.uniform(size=n), scale=1/shape)

    # Exponential to Beta: gbeta = exp(-gexp) ~ beta(shape, 1)
    gbeta = np.exp(-gexp)

    # Beta to Gamma: for V ~ gamma(1 + shape, rate),
    # X = gbeta*V ~ gamma(shape, rate) and Y = (1-gbeta)*V ~ gamma(1, rate)
    lV = np.log(gamma.rvs(1 + shape, scale=1/rate, size=n))
    Y = lV + np.vectorize(lambda j: logminus(1, j))(-gexp)
    X = lV - gexp

    return np.concatenate([X, Y])

# Function to compute qgamma_transformed
def qgamma_transformed(u1, u2, shape, rate): # pass
    # Uniform to Exponential
    gexp = expon.ppf(u1, scale=1/shape)

    # Exponential to Beta
    # gbeta = np.exp(-gexp)

    # Beta to Gamma
    V = gamma.ppf(u2, 1 + shape, scale=1/rate)
    lV = np.log(V)
    Y = lV + logminus(1, -gexp)
    X = lV - gexp

    return np.array([X, Y])


################
# transformed functions for dynesty
#################
def int_prior_transform(u): # pass
    # must have some initial points that u[xi] > 0.999
    if np.any(u < np.finfo(float).eps) or np.any(u >= 1):
        return -np.inf
    else:
        A = 2.5e07
        Time = 50000
        muloc = 1.2 * 10**(-6)
        thetaloc = 8.7 * 10**(-12)
        lmu0 = np.log(10 ** 6) - np.log(A) - np.log(Time)
        ltheta0 = np.log(10 ** 18) - 2 * (np.log(A) + np.log(Time))

        umu = np.array(LaplacesDemon.qtrunc(float(u[0]), spec = "cauchy", location = float(muloc), scale = float(muloc * 100), a = 0))
        utheta = np.array(LaplacesDemon.qtrunc(float(u[1]), spec = "cauchy", location = float(thetaloc), scale = float(thetaloc * 100), a = 0))
        upid = u[2]
        uxi = qgamma_bounded(u[3], shape=np.exp(2 * lmu0 - ltheta0), rate=np.exp(lmu0 - ltheta0))

        return np.concatenate((umu, utheta, np.array([upid, uxi])), axis = None)

def prior_transform(u): # pass
    # must have some initial points that u[xi] > 0.999
    if np.any(u < np.finfo(float).eps) or np.any(u >= 1):
        return -np.inf
    else:
        A = 2.5e07
        Time = 50000
        muloc = 1.2 * 10**(-6)
        thetaloc = 8.7 * 10**(-12)
        lmu0 = np.log(10**6) - np.log(A) - np.log(Time)
        ltheta0 = np.log(10**18) - 2 * (np.log(A) + np.log(Time))

        umu = np.array(LaplacesDemon.qtrunc(float(u[0]), spec = "cauchy", location = float(muloc), scale = float(muloc * 100), a = 0))
        utheta = np.array(LaplacesDemon.qtrunc(float(u[1]), spec = "cauchy", location = float(thetaloc), scale = float(thetaloc * 100), a = 0))
        upid = u[2]
        uxi = qgamma_bounded(u[3], shape=np.exp(2 * lmu0 - ltheta0), rate=np.exp(lmu0 - ltheta0))

        uLambda = q0inf_gamma(u[4:], shape=np.exp(2 * np.log(umu) - np.log(utheta)),
                               rate=np.exp(np.log(umu) - np.log(utheta)), pid=upid)

        return np.concatenate((umu, utheta, np.array([upid, uxi, 0]), uLambda), axis = None)

def log_dpois(loglambda, x): # pass
    if isinstance(x, (list, np.ndarray)):
        if isinstance(loglambda, (np.floating, np.integer, float, int)):
            return [x_i * loglambda - math.lgamma(x_i + 1) - np.exp(loglambda) for x_i in x]
        elif isinstance(loglambda, (list, np.ndarray)) and len(loglambda) == len(x):
            return [x_i * loglambda_i - math.lgamma(x_i + 1) - np.exp(loglambda_i) for x_i, loglambda_i in zip(x, loglambda)]
        else:
            raise ValueError("the lengths of loglambda and x do not match")
    elif isinstance(x, (np.floating, np.integer, float, int)):
        if isinstance(loglambda, (np.floating, np.integer, float, int)):
            return x * loglambda - math.lgamma(x + 1) - np.exp(loglambda)
        else:
            raise ValueError("Unsupported input types")
    else:
        raise ValueError("Unsupported input types")


def loglik_transform(t_xi, aux, Lambda, D, A, Time, a, r, e): # pass
    if t_xi >= -np.finfo(float).max and t_xi <= np.finfo(float).max:
        lfreq = np.log(A) + np.log(Time) + t_xi
        lfreqa = [Time * logplus(np.log(r) + np.log(e) + np.log(Lambda_i), np.log(a) + t_xi) for Lambda_i in Lambda]
        return log_dpois(loglambda=lfreq, x=D[-1]) + np.sum(log_dpois(loglambda=lfreqa, x=D[:-1]))
    else:
        return -1e100

def integrated_likelihood_transform(mu, theta, pid, t_xi, D, A = None, Time = None, a = None, R = None, e = None): # pass
    if A is None:
        A = simulation_functions.A
    if Time is None:
        Time = simulation_functions.Time
    if a is None:
        a = simulation_functions.a
    if R is None:
        R = simulation_functions.R
    if e is None:
        e = simulation_functions.e

    if mu <= 0 or theta <= 0 or pid < 0 or pid > 1 or t_xi <= -np.finfo(float).max or t_xi >= np.finfo(float).max:
        return -1e100
    else:
        t_mu = np.log(mu)
        t_theta = np.log(theta)
        t_pid = logit(pid)

        n = len(D) - 1
        shape1 = 1
        shape2 = 1
        muloc = 1.2e-6
        thetaloc = 8.7e-12
        mu0 = 10**6 / (A * Time)
        theta0 = 10**18 / (A * Time)**2
        Lamba_seq = np.arange(0, np.exp(7), 0.1)

        return analytical(t_mu=t_mu, t_theta=t_theta, t_xi=t_xi, t_pid=t_pid, D=D, A=A, Time=Time, a=a, r=R, e=e)


####################
# Posterior kernel
####################
import numpy as np
from scipy.stats import gamma, poisson, norm


def logpostker_gamma(mu, theta, pid, xi, Lambda, D, A, Time, a, r, e, mu0=None, theta0=None, muloc=None, thetaloc=None):
    # Default values if not provided
    if mu0 is None:
        mu0 = 1e6 / (A * Time)
    if theta0 is None:
        theta0 = 1e18 / (A * Time) ** 2
    if muloc is None:
        muloc = 1.2e-6
    if thetaloc is None:
        thetaloc = 8.7e-12

    return (loglik(xi, Lambda, D, A, Time, a, r, e) +
            logxiprior_gamma(xi, mu0, theta0) +
            logLambdaprior_gamma(Lambda, mu, theta, pid) +
            logpidhypprior(pid) +
            logmuhypprior(mu, muloc=muloc) +
            logthetahypprior(theta, thetaloc=thetaloc))


def t_logpostker_gamma(t_mu, t_theta, t_pid, t_xi, t_Lambda, D, A, Time, a, r, e,
                       mu0=None, theta0=None,
                       muloc=None, thetaloc=None):
    if mu0 is None:
        mu0 = 1e6 / (A * Time)
    if theta0 is None:
        theta0 = 1e18 / (A * Time) ** 2
    if muloc is None:
        muloc = 1.2e-6
    if thetaloc is None:
        thetaloc = 8.7e-12

    return (t_loglik(t_xi, t_Lambda, D, A, Time, a, r, e) +
            t_logxiprior_gamma(t_xi, mu0, theta0) +
            t_logLambdaprior_gamma(t_Lambda, t_mu, t_theta, t_pid) +
            t_logpidhypprior(t_pid) +
            t_logmuhypprior(t_mu, muloc=muloc) +
            t_logthetahypprior(t_theta, thetaloc=thetaloc))

def logpostker_unified_gamma(zeta, D, A, Time, a, R, e):
    def logpostker_in(zeta):
        return logpostker_gamma(mu=zeta[0],
                         theta=zeta[1],
                         pid=zeta[2],
                         xi=zeta[3],
                         Lambda=zeta[4:],
                         D=D, A=A, Time=Time, a=a, r=R, e=e)
    if np.ndim(zeta) == 1:
        return logpostker_in(zeta)
    elif np.ndim(zeta) == 2:
        return np.apply_along_axis(logpostker_in, axis=1, arr=zeta)

def t_logpostker_unified_gamma(zeta, D, A=None, Time=None, a=None, R=None, e=None):
    if A is None:
        A = simulation_functions.A
    if Time is None:
        Time = simulation_functions.Time
    if a is None:
        a = simulation_functions.a
    if R is None:
        R = simulation_functions.R
    if e is None:
        e = simulation_functions.e
    def t_logpostker_in(zeta):
        return t_logpostker_gamma(t_mu=zeta[0],
                           t_theta=zeta[1],
                           t_pid=zeta[2],
                           t_xi=zeta[3],
                           t_Lambda=zeta[4:],
                           D=D, A=A, Time=Time, a=a, r=R, e=e)

    if np.ndim(zeta) == 1:
        return t_logpostker_in(zeta)
    elif np.ndim(zeta) == 2:
        return np.apply_along_axis(t_logpostker_in, axis=1, arr=zeta)


