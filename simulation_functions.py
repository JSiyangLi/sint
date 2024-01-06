import numpy as np
from scipy.stats import truncexpon, gamma, beta, poisson, pareto
import LaplacesDemon

# Simulation 3.2 mh
# Static simulation settings
Time = 5e4
A = 2.5e7
xi = 2e-7
n = 10
break_point = 5


# Gamma lambdas
def Lambdastar_simulator(n, pid, mu_star=15, theta_star=50):
    # Step 2
    zero_indexes = np.where(np.log(np.random.uniform(size=n)) < np.log(pid))[0]
    Lambda_stars = gamma.rvs(size=n, a=mu_star ** 2 / theta_star, scale=mu_star / theta_star)
    Lambda_stars[zero_indexes] = 0
    return Lambda_stars


# Broken power law lambdas
def rbrokenpower(n, pllocation_vec, plshape_vec=None):
    if plshape_vec is None:
        plshape_vec = np.ones(len(pllocation_vec))

    def rmixtruncpareto_single(pllocation_vec, plshape_vec):
        pltau_diff = -plshape_vec[:-1] * np.diff(np.log(pllocation_vec))
        initial_logpjs = np.cumsum(np.concatenate(([np.logminus(x, 0) for x in pltau_diff], [0])))
        last_logpj = LaplacesDemon.logminus(0, LaplacesDemon.logplusvec(initial_logpjs))
        mix_weights = np.exp(np.concatenate(([initial_logpjs], [last_logpj]))) / \
                      np.sum(np.exp(np.concatenate(([initial_logpjs], [last_logpj]))))
        mix_truncparetos = np.array([truncexpon.rvs(size=1, loc=loc, scale=scale) for loc, scale in
                                     zip(pllocation_vec[:-1], pllocation_vec[1:])])
        mix_pareto = pareto.rvs(size=1, b=plshape_vec[-1], scale=pllocation_vec[-1])
        return np.sum(mix_weights * np.concatenate((mix_truncparetos, [mix_pareto])))

    return np.array([rmixtruncpareto_single(pllocation_vec, plshape_vec) for _ in range(n)])


def Lambdastar_simulator_brokenpower(pid, pllocation_vec, plshape_vec):
    # Step 2
    zero_indexes = np.where(np.log(np.random.uniform(size=n)) < np.log(pid))[0]
    Lambda_stars = rbrokenpower(n, pllocation_vec, plshape_vec)
    Lambda_stars[zero_indexes] = 0
    return Lambda_stars


# Simulator
def simulator(n, Lambda_stars, xi_star):
    # Step 1
    X = poisson.rvs(mu=A * xi * Time, size=1)
    # Step 3
    Background = poisson.rvs(mu=xi_star, size=n)
    Signal = poisson.rvs(mu=Lambda_stars, size=n)
    Y = Background + Signal
    return np.concatenate((Y, X))


# Prior sampling functions
def prior_sampling(sample_size=1, muloc=1.2e-6, thetaloc=8.7e-12, shape1=1, shape2=1, mu0=10 ** 6 / (A * Time),
                   theta0=10 ** 18 / (A * Time) ** 2):
    def qtcauchy(p, location=0, scale=1, a=None, b=None):
        second = 0 if a is None else beta.cdf(a, a, location=location, scale=scale)
        first = 1 if b is None else beta.cdf(b, a, location=location, scale=scale)
        return truncexpon.ppf(second + p * (first - second), location=location, scale=scale)

    def rtcauchy(n, location=0, scale=1, a=None, b=None):
        return qtcauchy(np.random.uniform(size=n), location, scale, a, b)

    mu_sample = rtcauchy(sample_size, location=muloc, scale=100 * muloc, a=0)
    theta_sample = rtcauchy(sample_size, location=thetaloc, scale=100 * thetaloc, a=0)
    pid_sample = beta.rvs(size=sample_size, a=shape1, b=shape2)

    uniflow = 1
    while uniflow >= 1e-300:
        tgamma = gamma.rvs(size=min(int(1e6), int(sample_size * 1e3)), a=mu0 ** 2 / theta0, scale=mu0 / theta0)
        uniflow = np.maximum(np.min(tgamma[tgamma != 0]), np.finfo(float).tiny)

    plowbound = gamma.cdf(uniflow, a=mu0 ** 2 / theta0, scale=mu0 / theta0, loc=0, scale=1, log=True)
    xi_sample = gamma.ppf(np.log(np.random.uniform(size=sample_size, low=uniflow, high=1)),
                          a=mu0 ** 2 / theta0, scale=mu0 / theta0, loc=0, scale=1)

    return np.column_stack((mu_sample, theta_sample, pid_sample, xi_sample))


# Test
prior_sampling(sample_size=10)

# Dynamic simulation settings
# Initial and true values
pid = 0.5
theta_star = 50
xi_star = 15
mu_star = 15
a = np.repeat(xi_star * 100, n)  # 100 from true values of Time and A
R = 1  # Modify for estimating Lambda, mu, theta
e = 1  # Modify for estimating Lambda, mu, theta
mu = mu_star / (R * e * Time)  # True value
theta = theta_star / (R * e * Time) ** 2  # True value
shape = mu ** 2 / theta
rate = mu / theta

# Transform initial values
t_mu = np.log(mu)
t_theta = np.log(theta)
t_pid = LaplacesDemon.logit(pid)
t_xi = np.log(xi)
t_shape = np.log(shape)
t_rate = np.log(rate)