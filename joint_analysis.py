import os
import numpy as np
import statsmodels.api as sm
import matplotlib
import rpy2.robjects.lib.ggplot2 as gp
import statsmodels.graphics.gofplots as sgplot
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
my_path = os.path.abspath(__file__)

emcee_flat_sample = np.load("full_emcee.npy")
emcee_int_flat_sample = np.load("marginal_emcee.npy")


from rpy2.robjects.packages import importr
utils = importr('utils')
base = importr('base')
stats = importr('stats')


#def foo(l, dtype=float):
#    return base.unlist(list((map(dtype, l))))

#param_name = ["ln(mu)", "ln(theta)", "logit(pid)", "ln(xi)"]
#print(list(foo(emcee_flat_sample[:, 0])), type(foo(emcee_flat_sample[:, 0])))
#p = (gp.ggplot(mapping=gp.aes(x = base.sort(foo(emcee_flat_sample[:, 0])), y = base.sort(foo(emcee_int_flat_sample[:, 0])))) +
#     gp.geom_point() +
#     gp.geom_abline(gp.aes(slope = 1, intercept = 0), linetype = 2))

def save_image(filename):
    # PdfPages is a wrapper around pdf
    # file so there is no clash and create
    # files with no error.
    p = PdfPages(filename)

    # get_fignums Return list of existing
    # figure numbers
    fig_nums = plt.get_fignums()
    figs = [plt.figure(n) for n in fig_nums]

    # iterating over the numbers in list
    for fig in figs:
        # and saving the files
        fig.savefig(p, format='pdf')

        # close the object
    p.close()

#emcee_qq = PdfPages("emcee_int_emcee")
for i in [0, 1, 2, 3]:
    pp_y = sm.ProbPlot(emcee_flat_sample[:, i], fit=True)
    pp_x = sm.ProbPlot(emcee_int_flat_sample[:, i], fit=True)
    fig = pp_x.qqplot(line="45", other=pp_y, xlabel="marginal emcee", ylabel = "full emcee")
    fig_nums = plt.get_fignums()
    figs = [plt.figure(n) for n in fig_nums]
    plt.show()
#emcee_qq.close()
#plt.savefig(os.path.join(my_path.replace("joint_analysis.py", ""), "py_plots/qq/emcee_int_emcee.pdf"), format = "pdf")


#sgplot.qqplot_2samples(emcee_flat_sample[:, 0], emcee_int_flat_sample[:, 0])
#plt.show()