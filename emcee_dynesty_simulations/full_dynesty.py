import numpy as np
from dynesty.utils import resample_equal
from dynesty import DynamicNestedSampler
import time
import sys
from kernel_functions import prior_transform, loglik
from simulation_functions import D, A, Time, a, R, e, n

# prior

def loglikelihood_dynesty(zeta):
    #mu = zeta[0]
    #theta = zeta[1]
    #pid = zeta[2]
    t_xi = zeta[3]
    #aux = zeta[4]
    Lambda = zeta[4:]
    output = loglik(t_xi, Lambda, D, A, float(Time), a, R, e)
    return output

nlive = 1024      # number of (initial) live points
bound = 'multi'   # use MutliNest algorithm
ndims = int(n + 4)
np.random.seed(42)

dsampler = DynamicNestedSampler(loglikelihood_dynesty, prior_transform, ndims,
                                    bound=bound)
start = time.time()
dsampler.run_nested(nlive_init=nlive)
dynesty_run_time = time.gmtime(time.time() - start)

dres = dsampler.results
dynesty_sample = dres.samples


dlogZdynesty = dres.logz[-1]        # value of logZ
dlogZerrdynesty = dres.logzerr[-1]  # estimate of the statistcal uncertainty on logZ

# output marginal likelihood
print('Marginalised evidence (using dynamic sampler) is {} ± {}'.format(dlogZdynesty, dlogZerrdynesty))

# get the posterior samples
dweights = np.exp(dres['logwt'] - dres['logz'][-1])
dpostsamples = resample_equal(dres.samples, dweights)

print('Number of posterior samples (using dynamic sampler) is {}'.format(dpostsamples.shape[0]))

# plot posterior samples (if corner.py is installed)
try:
    import matplotlib as mpl
    mpl.use("Agg") # force Matplotlib backend to Agg
    import corner # import corner.py
except ImportError:
    sys.exit(1)

fig = corner.corner(dpostsamples, labels=[r"$m$", r"$c$"], truths=[m, c], hist_kwargs={'density': True})
fig = corner.corner(dynesty_sample, fig=fig, color='r', hist_kwargs={'density': True})
