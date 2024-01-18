import numpy as np
from simulation_functions import D, t_mu, t_theta, t_pid, t_xi, truths, labels
from kernel_functions import t_int_kernel_gamma
import emcee
import time
from tqdm import tqdm
import h5py
import matplotlib as mpl
mpl.interactive(True)
import matplotlib.pyplot as plt
import corner # import corner.py
import os
my_path = os.path.abspath(__file__)

#############################
# emcee parameters
#######################

# sampler settings
np.random.seed(40)
nwalkers = int(8)
ndim = 4
emcee_int_sampler = emcee.EnsembleSampler(
    nwalkers=nwalkers,
    ndim=ndim,
    log_prob_fn=t_int_kernel_gamma,
    vectorize=True,
    args = [D]
)

##################
# initial values
#################
turb_scale = 2
turb_mu = np.random.normal(size = nwalkers, loc = t_mu, scale = np.abs((t_mu + 1) / turb_scale))
turb_mu = np.expand_dims(np.array(turb_mu), axis = -1)
turb_theta = np.random.normal(size = nwalkers, loc = t_theta, scale = np.abs((t_theta + 1) / turb_scale))
turb_theta = np.expand_dims(np.array(turb_theta), axis = -1)
turb_pid = np.random.normal(size = nwalkers, loc = t_pid, scale = np.abs((t_pid + 1) / turb_scale))
turb_pid = np.expand_dims(np.array(turb_pid), axis = -1)
turb_xi = np.random.normal(size = nwalkers, loc = t_xi, scale = np.abs((t_xi + 1) / turb_scale))
turb_xi = np.expand_dims(np.array(turb_xi), axis = -1)
t_int_init = np.hstack((turb_mu, turb_theta, turb_pid, turb_xi))

# testing the posterior kernel function in python
t_int_kernel_gamma(zeta=t_int_init, D=D)
t_int_kernel_gamma(t_int_init[[1], :], D=D)

# warnings.filterwarnings("ignore")
# HDFbackend = emcee.backends.HDFBackend("emcee_chains")
start = time.time()
emcee_int_result = emcee_int_sampler.run_mcmc(t_int_init, nsteps=50000, progress=True,
                                              skip_initial_state_check=True)
int_emcee_run_time = time.time() - start
print(int_emcee_run_time)

emcee_int_sample = emcee_int_sampler.get_chain()
t_int_kernel_gamma(zeta=emcee_int_result.coords, D=D)
thin = 10
emcee_int_flat_sample = emcee_int_sampler.get_chain(flat=True, thin = thin)
int_emcee_autoco = emcee_int_sampler.get_autocorr_time()

np.mean(emcee_int_flat_sample[:, 1])
fig = corner.corner(emcee_int_flat_sample, labels=labels, truths=truths, hist_kwargs={'density': True})
plt.savefig(os.path.join(my_path.replace("marginal_emcee.py", ""), "py_plots/corner/marginal_emcee.pdf"), format = "pdf")

######################
# save the posterior
#####################
np.save("marginal_emcee.npy", emcee_int_flat_sample)