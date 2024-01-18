import numpy as np
from dynesty.utils import resample_equal
from dynesty import NestedSampler, DynamicNestedSampler
from kernel_functions import int_prior_transform, integrated_likelihood_transform
from simulation_functions import D, A, Time, a, R, e, truths, labels
import time, sys
import matplotlib as mpl
from matplotlib import pyplot as plt
mpl.use("Agg") # force Matplotlib backend to Agg
import corner # import corner.py
from dynesty import plotting as dyplot
from rpy2.robjects.packages import importr
utils = importr('utils')
base = importr('base')
stats = importr('stats')
LaplacesDemon = importr('LaplacesDemon')

np.random.seed(42)
# prior
def int_loglikelihood_dynesty(zeta):
    mu = zeta[0]
    theta = zeta[1]
    pid = zeta[2]
    xi = zeta[3]
    return integrated_likelihood_transform(mu, theta, pid, xi, D, A, float(Time), a, R, e)

int_nlive = 1024      # number of (initial) live points
int_bound = 'multi'   # use MutliNest algorithm
int_ndims = int(4)

int_dsampler = DynamicNestedSampler(int_loglikelihood_dynesty, int_prior_transform, int_ndims,
                                    bound=int_bound)
start = time.time()
int_dsampler.run_nested(nlive_init=int_nlive)
int_dynesty_run_time = time.gmtime(time.time() - start)
int_dres = int_dsampler.results
int_dynesty_sample = int_dres.samples


dlogZdynesty = int_dres.logz[-1]        # value of logZ
dlogZerrdynesty = int_dres.logzerr[-1]  # estimate of the statistcal uncertainty on logZ

# output marginal likelihood
print('Marginalised evidence (using dynamic sampler) is {} ± {}'.format(dlogZdynesty, dlogZerrdynesty))

# get the posterior samples
int_dweights = np.exp(int_dres['logwt'] - int_dres['logz'][-1])
int_dpostsamples = resample_equal(int_dres.samples, int_dweights)
t_int_dpostsamples = [np.log(int_dpostsamples.samples[:, 0]), np.log(int_dpostsamples.samples[:, 1]), np.array(LaplacesDemon.logit(int_dpostsamples.samples[:, 1])), np.log(int_dpostsamples.samples[:, 3])]

print('Number of posterior samples (using dynamic sampler) is {}'.format(int_dpostsamples.shape[0]))

# plot posterior samples (if corner.py is installed)
fig = corner.corner(t_int_dpostsamples, labels=labels, truths=truths, hist_kwargs={'density': True})
plt.show()
