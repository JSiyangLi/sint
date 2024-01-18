from kernel_functions import t_logpostker_unified_gamma
from simulation_functions import n, mu, theta, pid, xi, t_mu, t_theta, t_pid, t_xi, D, labels, truths

import emcee
import time
import numpy as np
from tqdm import tqdm
import h5py
import matplotlib as mpl
mpl.interactive(True)
import matplotlib.pyplot as plt
import corner # import corner.py
import os
my_path = os.path.abspath(__file__)


from rpy2.robjects.packages import importr
utils = importr('utils')
base = importr('base')
stats = importr('stats')

# sampler settings
np.random.seed(50)
nwalkers = int(2 * (n + 4))
ndim = n + 4

emcee_sampler = emcee.EnsembleSampler(
  nwalkers = nwalkers,
  ndim = ndim,
  log_prob_fn = t_logpostker_unified_gamma,
  vectorize = True,
  args = [D]
)

# original scale true values
truevals = np.array([mu, theta, pid, xi] * nwalkers)
truevals = np.multiply(truevals.reshape((int(nwalkers), 4)), np.random.rand(nwalkers, 4))
lambs = np.random.gamma(shape = (mu**2) / theta, scale = theta / mu, size = n * nwalkers)
lambs = np.expand_dims(np.array(lambs), axis = -1)
lambs = lambs.reshape((int(nwalkers), int(n)))
init = np.hstack((truevals, lambs))

#initmean = np.mean(init, axis=0)
#initstd = np.std(init, axis=0)
#init_standardized = (init - initmean) / initstd

# transformed scale true values
# hyperparameters
turb_scale = 2
turb_mu = np.random.normal(size = nwalkers, loc = t_mu, scale = np.abs((t_mu + 1) / turb_scale))
turb_mu = np.expand_dims(np.array(turb_mu), axis = -1)
turb_theta = np.random.normal(size = nwalkers, loc = t_theta, scale = np.abs((t_theta + 1) / turb_scale))
turb_theta = np.expand_dims(np.array(turb_theta), axis = -1)
turb_pid = np.random.normal(size = nwalkers, loc = t_pid, scale = np.abs((t_pid + 1) / turb_scale))
turb_pid = np.expand_dims(np.array(turb_pid), axis = -1)
turb_xi = np.random.normal(size = nwalkers, loc = t_xi, scale = np.abs((t_xi + 1) / turb_scale))
turb_xi = np.expand_dims(np.array(turb_xi), axis = -1)

# lambdas
#large_random_order = 1.1
# for the ones with lambda > 0
#t_lambs_indicator = Lambda_stars == 0
#t_lambs = r.rnorm(n = nwalkers * r.n, mean = np.max(Lambda_stars), sd = np.abs(np.max(Lambda_stars + 1) / turb_scale))
#t_lambs = np.expand_dims(np.array(t_lambs), axis = -1)
#t_lambs = t_lambs.reshape((int(nwalkers), int(r.n)))
# for the ones with lambda = 0 (true value = -infinity)
#large_neg_random = -np.random.rand(int(nwalkers * np.sum(t_lambs_indicator))) * 10**(large_random_order * np.random.rand(int(nwalkers * np.sum(t_lambs_indicator)))) # order should not exceed the order of sys.float_info.max
#large_neg_random = np.expand_dims(np.array(large_neg_random), axis = -1)
#large_neg_random = large_neg_random.reshape((int(nwalkers), np.sum(t_lambs_indicator)))
#t_lambs[:, t_lambs_indicator] = large_neg_random

t_lambs = np.log(np.random.gamma(size = n * nwalkers, shape = mu**2 / theta, scale = theta / mu))
t_lambs = np.expand_dims(np.array(t_lambs), axis = -1)
t_lambs = t_lambs.reshape((int(nwalkers), int(n)))

t_init = np.hstack((turb_mu, turb_theta, turb_pid, turb_xi, t_lambs))

# testing the posterior kernel function in python
#print("zeta", t_init)

#print(t_logpostker_unified_gamma(zeta = t_init, D = D))
print(t_logpostker_unified_gamma(t_init[[1], :], D = D))


#warnings.filterwarnings("ignore")
#HDFbackend = emcee.backends.HDFBackend("emcee_chains")
start = time.time()
emcee_result = emcee_sampler.run_mcmc(t_init, nsteps = 50000, progress = True,
  skip_initial_state_check = True)
emcee_run_time = time.gmtime(time.time() - start)
print(emcee_run_time)

emcee_sample = emcee_sampler.get_chain()
#t_logpostker_unified_gamma(zeta = emcee_result.coords, D = D)
thin = 10  # Adjust as needed
emcee_flat_sample = emcee_sampler.get_chain(flat=True, thin=thin)
emcee_autoco = emcee_sampler.get_autocorr_time()

fig = corner.corner(emcee_flat_sample[:, [0, 1, 2, 3]], labels=labels, truths=truths, hist_kwargs={'density': True})
plt.savefig(os.path.join(my_path.replace("full_emcee.py", ""), "py_plots/corner/full_emcee.pdf"), format = "pdf")

######################
# save the posterior
#####################
np.save("full_emcee.npy", emcee_flat_sample)