import numpy as np
from scipy.stats import cauchy, beta, gamma
#######################
# simulation 3.2 mh
###########
# static simulation settings
Time = 5e4
A = 2.5e7
xi = 2e-7
n = 10
break_point = 5

def Lambdastar_simulator(n, pid, mu_star=15, theta_star=1):
    # Step 2
    zero_indexes = np.where(np.log(np.random.uniform(size=n)) < np.log(pid))
    Lambda_stars = np.random.gamma(shape=mu_star ** 2 / theta_star, scale=theta_star / mu_star, size=n)
    Lambda_stars[zero_indexes] = 0

    return Lambda_stars


def simulator(n, Lambda_stars, xi_star, A, Time):
    # Step 1
    X = np.random.poisson(A * xi * Time, size=1)

    # Step 3
    Background = np.random.poisson(xi_star, size=n)
    Signal = np.random.poisson(Lambda_stars, size=n)
    Y = Background + Signal

    return np.concatenate((Y, X))

def prior_sampling(sample_size=1, muloc=1.2e-6, thetaloc=8.7e-12, shape1=1, shape2=1, mu0=10 ** 6, theta0=10 ** 18):
    # Truncated Cauchy
    def qtcauchy(p, location=0, scale=1, a=None, b=None):
        second = 0 if a is None else cauchy.cdf(a, loc=location, scale=scale)
        first = 1 if b is None else cauchy.cdf(b, loc=location, scale=scale)
        return cauchy.ppf(second + p * (first - second), loc=location, scale=scale)

    def rtcauchy(n, location=0, scale=1, a=None, b=None):
        return qtcauchy(np.random.uniform(size=n), location, scale, a, b)

    # mu
    mu_sample = rtcauchy(sample_size, location=muloc, scale=100 * muloc, a=0)

    # theta
    theta_sample = rtcauchy(sample_size, location=thetaloc, scale=100 * thetaloc, a=0)

    # pid
    pid_sample = beta.rvs(shape1, shape2, size=sample_size)

    # xi
    uniflow = 1
    while uniflow >= 1e-300:
        tgamma = gamma.rvs((mu0 ** 2) / theta0, loc=0, scale=mu0 / theta0, size=min(int(1e6), int(sample_size * 1e3)))
        uniflow = max(np.min(tgamma[tgamma != 0]), np.finfo(float).tiny)  # Avoiding division by zero
    plowbound = gamma.logcdf(uniflow, (mu0 ** 2) / theta0, loc=0, scale=mu0 / theta0)
    xi_sample = gamma.ppf(np.exp(np.log(np.random.uniform(size=sample_size, low=np.exp(plowbound), high=1))),
                          (mu0 ** 2) / theta0, loc=0, scale=mu0 / theta0)

    return np.column_stack((mu_sample, theta_sample, pid_sample, xi_sample))

##############################
# dynamic simulation settings
##############################
# gamma
######################################
# initial and true values
pid = 0.5
theta_star = 50
xi_star = 15
mu_star = 15
a = xi_star * 100 # 100 from true values of Time and A
R = 1 # modify for estimating Lambda, mu, theta
e = 1 # modify for estimating Lambda, mu, theta
mu = mu_star / (R * e * Time) # true value
theta = theta_star / (R * e * Time)**2 # true value
shape = mu**2 / theta
rate = mu / theta

# transform initial values
t_mu = np.log(mu)
t_theta = np.log(theta)
t_pid = np.log(pid / (1 - pid))
t_xi = np.log(xi)
t_shape = np.log(shape)
t_rate = np.log(rate)
truths = [t_mu, t_theta, t_pid, t_xi]
labels = [r"$ln(mu)$", r"$ln(theta)$", r"$logit(pid)$", r"$ln(xi)$"]

##############
# simulator
############
np.random.seed(0)
Lambda_stars = np.array(Lambdastar_simulator(n = n, pid = pid, theta_star = theta_star))
D = np.array(simulator(n = n, Lambda_stars = Lambda_stars, xi_star = xi_star, A = A, Time = Time))