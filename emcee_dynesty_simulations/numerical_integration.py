from scipy.stats import poisson, gamma
from scipy.special import gammaln, expit
import sys
import numpy as np
from rpy2.robjects.packages import importr
utils = importr('utils')
base = importr('base')
stats = importr('stats')

def logplus(x, y):
    if x > y:
        return x + np.log(1.0 + np.exp(y - x))
    else:
        return y + np.log(1.0 + np.exp(x - y))

def logplusvec(x):
    n = len(x)
    r = -sys.float_info.max

    for i in range(n):
        r = logplus(r, x[i])

    return r

def logminus(x, y):
    if x >= y:
        return x + np.log(1.0 - np.exp(y - x))
    else:
        return float('nan')

def logdeltafunction(i): # pass
    return 0 if i == 0 else -np.inf


def t_Xlik(x, t_xi, A, Time): # pass
    xi = np.exp(t_xi)
    return poisson.logpmf(x, A * Time * xi)


def t_Lambda_ker(y, Lambda, t_mu, t_theta, t_pid, t_xi, A, Time, a, r, e): # pass
    pid = 1 / (1 + np.exp(-t_pid))
    xi = np.exp(t_xi)
    # Jacobian calculation is not needed for Lambda
    def prior_func(y):
        return logplus(np.log(pid) + logdeltafunction(y),
                       np.log(1 - pid) + gamma.logpdf(y, a=np.exp(2 * t_mu - t_theta), scale=np.exp(t_theta - t_mu), loc=0))
    lprior = [prior_func(j) for j in Lambda]
    llike = [poisson.logpmf(y_l, (a * xi + r * e * Lambda_l) * Time) for y_l, Lambda_l in zip(y, Lambda)]
    output = np.add(lprior, llike)
    return output


def t_Lambda_gamma(y, Lambda, t_mu, t_theta, t_xi, Time, a, r, e): # pass
    lshape = 2 * t_mu - t_theta
    lrate = t_mu - t_theta
    xi = np.exp(t_xi)

    if not isinstance(Lambda, list):
        Poisrate = float((a * xi + r * e * Lambda) * Time)
    else:
        Poisrate = [(a * xi + r * e * Lambda_j) * Time_j for Lambda_j, Time_j in zip(Lambda, Time)]

    lprior = gamma.logpdf(Lambda, a=np.exp(lshape), scale=np.exp(-lrate), loc=0)
    llike = stats.dpois(float(y), Poisrate, log = True)
    return np.add(lprior, np.array(llike))

#########
# numerical integration
#########
def numerical(t_mu, t_theta, t_pid, t_xi, Lambda, x, y, A, Time, a, r, e): # pass
    pid = 1 / (1 + np.exp(-t_pid))
    xi = np.exp(t_xi)

    def ut_Lambda_gamma(y, Lambda):
        return t_Lambda_gamma(y, Lambda, t_mu, t_theta, t_xi, Time, a, r, e)
    u_Lambda_gamma = np.frompyfunc(ut_Lambda_gamma, 2, 1)
    o = u_Lambda_gamma.outer(y, [num for num in Lambda if num]) # the output is a transpose of the outer() by R
    into = np.apply_along_axis(logplusvec, axis=1, arr=o) - np.log(o.shape[1])

    zeroinf = np.log(pid) + poisson.logpmf(y, a * xi * Time)
    nonzeroinf = np.log(1 - pid) + into

    return np.add(t_Xlik(x, t_xi, A, Time), np.sum([logplus(zero_i, nonzero_i) for zero_i, nonzero_i in zip(zeroinf, nonzeroinf)]))

def analytical(t_mu=None, t_theta=None, t_pid=None, t_xi=None,
               D=None, A=None, Time=None, a=None, r=None, e=None): # pass
    # reshape formal arguments
    X = D[-1]
    Y = D[:-1]

    t_shape = 2 * t_mu - t_theta
    t_rate = t_mu - t_theta
    # construction
    def summand(yi, k, t_shape, t_rate, t_xi, Time, a, r, e):
        lcombination = gammaln(yi + 1) - (gammaln(k + 1) + gammaln(yi - k + 1))
        lpowers = (yi - k) * (np.log(a) + t_xi + np.log(Time)) + k * (np.log(r) + np.log(e) + np.log(Time))
        llaplace = gammaln(np.exp(t_shape) + k) + (-np.exp(t_shape) - k) * (np.logaddexp(t_rate, np.log(r) + np.log(e) + np.log(Time)))
        return lcombination + lpowers + np.where((np.exp(t_shape) + k) > 0, llaplace, -np.inf)  # the condition of Laplace transform

    def prodant(yi, t_shape, t_rate, t_xi, t_pid, Time, a, r, e):
        pid = expit(t_pid)

        lpois = -a * np.exp(t_xi) * Time - gammaln(yi + 1)
        l0inf = yi * (np.log(a) + t_xi + np.log(Time)) + np.log(pid)
        sums = np.log(1 - pid) + np.exp(t_shape) * t_rate - gammaln(np.exp(t_shape)) + \
            logplusvec(summand(yi, np.arange(yi + 1), t_shape, t_rate, t_xi, Time, a, r, e))
        return lpois + logplus(l0inf, sums)

    def prodant_vec(Y, t_shape, t_rate, t_xi, t_pid, Time, a, r, e):
        return np.array([prodant(yi, t_shape, t_rate, t_xi, t_pid, Time, a, r, e) for yi in Y])

    return (np.sum(prodant_vec(Y, t_shape, t_rate, t_xi, t_pid, Time, a, r, e)) +
            poisson.logpmf(X, A * Time * np.exp(t_xi)))

