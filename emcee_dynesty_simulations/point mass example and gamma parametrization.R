library(Rcpp)
library(matrixStats)
sourceCpp("logsum.cpp")
#####################
# example 1: Unif-normal - large point mass prob
set.seed(42)
n = 1000
point_mass_prob <- 0.9
X_max = 4
X <- runif(n = n, max = X_max)
replacement <- runif(n = n) < point_mass_prob # controls how many points are dominated by the normal component.
X[replacement] <- 0
Y <- X + rnorm(n = n)
plot(density(Y),
     main = "Density of rnorm based on Unif-point-mass mixture")
# numerical integral over X given Y values:
integrand_interval <- seq(0, 4, by = 0.01)
individual_obs <- point_mass_prob * dnorm(Y, mean = 0) + # integrate step:
  (1 - point_mass_prob) * colSums(sapply(Y, function(y) 0.25 * dnorm(y, mean = integrand_interval))) / length(integrand_interval)
tot_mix <- sum(log(individual_obs))
hist(individual_obs, main = "histogram of integrand for each observation")
pdf("meeting_shared/9Nov2023/Uniform example.pdf")
plot(density(individual_obs), main = "density of log integrand for each observation",
     sub = paste0("point mass prob = ", point_mass_prob, ", total mixture = ", tot_mix))

# small point mass prob
set.seed(42)
n = 1000
point_mass_prob <- 0.1
X_max = 4
X <- runif(n = n, max = X_max)
replacement <- runif(n = n) < point_mass_prob # controls how many points are dominated by the normal component.
X[replacement] <- 0
Y <- X + rnorm(n = n)

# numerical integral over X given Y values:
integrand_interval <- seq(0, 4, by = 0.01)
individual_obs <- point_mass_prob * dnorm(Y, mean = 0) + # integrate step:
  (1 - point_mass_prob) * colSums(sapply(Y, function(y) 0.25 * dnorm(y, mean = integrand_interval))) / length(integrand_interval)
tot_mix <- sum(log(individual_obs))
#hist(individual_obs, main = "histogram of integrand for each observation")
plot(density(individual_obs), main = "density of integrand for each observation",
     sub = paste0("point mass prob = ", point_mass_prob, ", total mixture = ", tot_mix))
dev.off()

###############################
# example 2: bivariate normal - large point mass prob, rho = 0.9
#library(ks)
pdf("meeting_shared/9Nov2023/Normal example.pdf")
set.seed(42)
n = 1000
rho <- 0.9
point_mass_prob = 0.9
X <- rnorm(n = n)
replacement <- runif(n = n) < point_mass_prob # controls how many points are dominated by the point mass on x.
X[replacement] <- 2
Y <- rnorm(n = n, mean = 3 * X * rho, sd = sqrt(1 - rho^2))
#plot(density(Y))
plot(X, Y, pch = 19, main = paste0("Point mass prob = ", point_mass_prob, ", rho = ", rho))

#cont = seq(5, 95, by = 10)
#XY_kde <- kde(cbind(X, Y))
#plot(XY_kde, display = "filled.contour", cont = cont, col.fun = viridis::viridis,
#     xlab = "X", ylab = "Y")
#plot(XY_kde, cont = cont, display = "slice", add = TRUE)

# numerical integral over X given Y values:
integrand_interval <- seq(-15, 15, by = 0.01)
individual_obs <- point_mass_prob * dnorm(Y, mean = 2) + # integrate step:
  (1 - point_mass_prob) * colSums(sapply(Y, function(y) dnorm(integrand_interval) * dnorm(y, mean = integrand_interval))) / length(integrand_interval)
tot_mix <- sum(log(individual_obs))
#hist(individual_obs, main = "histogram of integrand for each observation")
plot(density(individual_obs), main = "Density of integrand for each observation",
     sub = paste0("point mass prob = ", point_mass_prob, ", total mixture = ", tot_mix))
mtext(text = paste0("rho = ", rho))

# small point mass prob, rho = 0.9
set.seed(42)
n = 1000
rho <- 0.9
point_mass_prob = 0.1
X <- rnorm(n = n)
replacement <- runif(n = n) < point_mass_prob # controls how many points are dominated by the point mass on x.
X[replacement] <- 2
Y <- rnorm(n = n, mean = 3 * X * rho, sd = sqrt(1 - rho^2))
#plot(density(Y))
plot(X, Y, pch = 19, main = paste0("Point mass prob = ", point_mass_prob, ", rho = ", rho))

#cont = seq(5, 95, by = 10)
#XY_kde <- kde(cbind(X, Y))
#plot(XY_kde, display = "filled.contour", cont = cont, col.fun = viridis::viridis,
#     xlab = "X", ylab = "Y")
#plot(XY_kde, cont = cont, display = "slice", add = TRUE)

# numerical integral over X given Y values:
integrand_interval <- seq(-15, 15, by = 0.01)
individual_obs <- point_mass_prob * dnorm(Y, mean = 2) + # integrate step:
  (1 - point_mass_prob) * colSums(sapply(Y, function(y) dnorm(integrand_interval) * dnorm(y, mean = integrand_interval))) / length(integrand_interval)
tot_mix <- sum(log(individual_obs))
#hist(individual_obs, main = "histogram of integrand for each observation")
plot(density(individual_obs), main = "Density of integrand for each observation",
     sub = paste0("point mass prob = ", point_mass_prob, ", total mixture = ", tot_mix))
mtext(text = paste0("rho = ", rho))

# large point mass prob, rho = 0.1
set.seed(42)
n = 1000
rho <- 0.1
point_mass_prob = 0.9
X <- rnorm(n = n)
replacement <- runif(n = n) < point_mass_prob # controls how many points are dominated by the point mass on x.
X[replacement] <- 2
Y <- rnorm(n = n, mean = 3 * X * rho, sd = sqrt(1 - rho^2))
#plot(density(Y))
plot(X, Y, pch = 19, main = paste0("Point mass prob = ", point_mass_prob, ", rho = ", rho))

#cont = seq(5, 95, by = 10)
#XY_kde <- kde(cbind(X, Y))
#plot(XY_kde, display = "filled.contour", cont = cont, col.fun = viridis::viridis,
#     xlab = "X", ylab = "Y")
#plot(XY_kde, cont = cont, display = "slice", add = TRUE)

# numerical integral over X given Y values:
integrand_interval <- seq(-15, 15, by = 0.01)
individual_obs <- point_mass_prob * dnorm(Y, mean = 2) + # integrate step:
  (1 - point_mass_prob) * colSums(sapply(Y, function(y) dnorm(integrand_interval) * dnorm(y, mean = integrand_interval))) / length(integrand_interval)
tot_mix <- sum(log(individual_obs))
#hist(individual_obs, main = "histogram of integrand for each observation")
plot(density(individual_obs), main = "Density of integrand for each observation",
     sub = paste0("point mass prob = ", point_mass_prob, ", total mixture = ", tot_mix))
mtext(text = paste0("rho = ", rho))

# small point mass prob, rho = 0.1
set.seed(42)
n = 1000
rho <- 0.1
point_mass_prob = 0.1
X <- rnorm(n = n)
replacement <- runif(n = n) < point_mass_prob # controls how many points are dominated by the point mass on x.
X[replacement] <- 2
Y <- rnorm(n = n, mean = 3 * X * rho, sd = sqrt(1 - rho^2))
#plot(density(Y))
plot(X, Y, pch = 19, main = paste0("Point mass prob = ", point_mass_prob, ", rho = ", rho))

#cont = seq(5, 95, by = 10)
#XY_kde <- kde(cbind(X, Y))
#plot(XY_kde, display = "filled.contour", cont = cont, col.fun = viridis::viridis,
#     xlab = "X", ylab = "Y")
#plot(XY_kde, cont = cont, display = "slice", add = TRUE)

# numerical integral over X given Y values:
integrand_interval <- seq(-15, 15, by = 0.01)
individual_obs <- point_mass_prob * dnorm(Y, mean = 2) + # integrate step:
  (1 - point_mass_prob) * colSums(sapply(Y, function(y) dnorm(integrand_interval) * dnorm(y, mean = integrand_interval))) / length(integrand_interval)
tot_mix <- sum(log(individual_obs))
#hist(individual_obs, main = "histogram of integrand for each observation")
plot(density(individual_obs), main = "Density of integrand for each observation",
     sub = paste0("point mass prob = ", point_mass_prob, ", total mixture = ", tot_mix))
mtext(text = paste0("rho = ", rho))
dev.off()

###########################
# example 3: 
source(paste0(getwd(), "/simulation_functions.R"))
source(paste0(getwd(), "/kernel_functions.R"))
emcee_results_list <- readRDS(file = paste0(getwd(), "/simple/emcee/results/emcee_result.RDS"))
# testing intervals and parameter true values
tmu_seq <- seq(-6, 0, by = 0.1)
ttheta_seq <- seq(-60, 0, by = 0.1)
tpid_seq <- seq(-10, 10, by = 0.01)
txi_seq <- seq(-18, -13, by = 0.01)
beta_seq <- alpha_seq <- seq(-50, 50, by = 0.1)
t_alpha <- log(4.5)
t_beta <- log(15000)

# on the mu-theta scale:
plot(tmu_seq,
     sapply(tmu_seq, function(k) integrated_kernel(t_mu = k, t_theta = t_theta, t_pid = t_pid, t_xi = t_xi, D = c(0:25, 30), A = A, Time = Time, a = a, r = R, e = e)),
     type = "l")
plot(ttheta_seq,
     sapply(ttheta_seq, function(k) integrated_kernel(t_mu = t_mu, t_theta = k, t_pid = t_pid, t_xi = t_xi, D = c(0:3, 30), A = A, Time = Time, a = a, r = R, e = e)),
     type = "l",
     xlab = "log theta", ylab = "conditional kernel", main = "Lazhi parametrized conditional Gamma kernel on theta") # not log concave. the dominating part is:
abline(v = t_theta, col = "red")
plot(ttheta_seq,
     sapply(ttheta_seq, function(k) {
       lshape <- 2 * t_mu - k
       exp(lshape) * (t_mu - k) - lgamma(exp(lshape)) + logplusvec(lgamma(exp(lshape) + 0:10) + (-exp(t_shape) - 0:10) * (logplus(t_rate, log(R) + log(e) + log(Time))))
     }),
     type = "l",
     xlab = "log theta", ylab = "Laplace(Gamma)", main = "Lazhi parametrized Gamma on theta")
abline(v = t_theta, col = "red")
# which is equivalent to:
plot(ttheta_seq, 
     sapply(ttheta_seq, function(k) {
       lshape <- 2 * t_mu - k
       lrate <- t_mu - k
       dgamma(1, shape = exp(lshape), rate = exp(lrate), log = TRUE) + 1 * exp(lrate) - (exp(lshape) - 1) * log(1)
     }),
     type = "l")

# on the shape-rate scale:
plot(beta_seq,
     exp(t_alpha) * beta_seq, #- 100 * exp(beta_seq),
     type = "l")
plot(beta_seq, 
     dgamma(100, shape = exp(t_alpha), rate = exp(beta_seq), log = TRUE),
     type = "l")
plot(alpha_seq,
     exp(alpha_seq) * t_beta - lgamma(exp(alpha_seq)), #+ (exp(alpha_seq) - 1) * log(100) - 100 * exp(t_beta),
     type = "l")
plot(alpha_seq, 
     dgamma(100, shape = exp(alpha_seq), rate = exp(t_beta), log = TRUE),
     type = "l")
# which is log concave, hence this solves the problem.

# after the reparametrization:
tshape_seq <- seq(-200, 200, by = 1)
trate_seq <- seq(-200, 200, by = 1)

# the rate parameter was fine
plot(trate_seq,
     sapply(trate_seq, function(k) integrated_kernel(t_shape, k, t_pid, t_xi, D = c(0:30, 30), A = A, Time = Time, a = a, r = R, e = e)),
     type = "l")
abline(v = t_rate, col = "green")
plot(trate_seq,
     sapply(trate_seq, function(k) analytical(t_shape, k, t_pid, t_xi, D = c(0:30, 30), A = A, Time = Time, a = a, r = R, e = e)),
     type = "l")
abline(v = t_rate, col = "green")
plot(trate_seq,
     sapply(trate_seq, function(k) t_logshaperatehypprior(t_shape, k)),
     type = "l") # this is how mu-theta-transformed Cauchy looks on the log scale, by nature.
abline(v = log(1.2e-6) - log(8.7e-12), col = "blue")
# the shape parameter is not concave enough
plot(tshape_seq,
     sapply(tshape_seq, function(k) integrated_kernel(k, t_rate, t_pid, t_xi, D = c(0:25, 30), A = A, Time = Time, a = a, r = R, e = e)),
     type = "l",
     xlab = "shape parameter alpha", ylab = "marginal posterior kernel", 
     main = "marginal posterior kernel vs. Gamma alpha")
abline(v = t_shape, col = "green")
plot(tshape_seq[-1],
     diff(sapply(tshape_seq, function(k) integrated_kernel(k, t_rate, t_pid, t_xi, D = c(0:3, 30), A = A, Time = Time, a = a, r = R, e = e))),
     type = "l")
abline(v = t_shape, col = "green")
plot(tshape_seq,
     sapply(tshape_seq, function(k) analytical(k, t_rate, t_pid, t_xi, D = c(0:30, 30), A = A, Time = Time, a = a, r = R, e = e)),
     type = "l",
     xlab = "shape parameter alpha", ylab = "integral", 
     main = "integral vs. Gamma alpha")
abline(v = t_shape, col = "green")
plot(tshape_seq,
     sapply(tshape_seq, function(p) summand(10, 0, p, t_rate, t_xi = t_xi, Time = Time, a = a, r = R, e = e)),
     type = "l",
     xlab = "shape parameter alpha", ylab = "main parts in integral", 
     main = "integral vs. Gamma alpha")
abline(v = t_shape, col = "green")
plot(tshape_seq,
     sapply(tshape_seq, function(k) t_logshaperatehypprior(k, t_rate)),
     type = "l") # this is how transformed Cauchy looks on the log scale, by nature.
abline(v = log(1.2e-6) - log(8.7e-12), col = "blue")

plot(ttheta_seq, peaky_avg_gamma[, 1], type = "l")

# saved plots
pdf("meeting_shared/9Nov2023/gamma example.pdf")
plot(ttheta_seq,
     sapply(ttheta_seq, function(k) integrated_kernel(t_mu = t_mu, 
                                                      t_theta = k,
                                                      t_pid = t_pid,
                                                      t_xi = t_xi, D = emcee_results_list$simulated_data, A = A, Time = Time, a = a, r = R, e = e)),
     type = "l",
     xlab = "log theta", ylab = "conditional kernel", main = "Lazhi parametrized conditional Gamma on theta") # not log concave. the dominating part is:
abline(v = t_theta, col = "red")
plot(ttheta_seq,
     sapply(ttheta_seq, function(k) integrated_kernel(t_mu = t_mu, 
                                                      t_theta = k,
                                                      t_pid = t_pid,
                                                      t_xi = t_xi, D = c(18, 30), A = A, Time = Time, a = a, r = R, e = e)),
     type = "l",
     xlab = "log theta", ylab = "conditional kernel", main = "Lazhi parametrized conditional Gamma on theta") # not log concave. the dominating part is:
abline(v = t_theta, col = "red")
lines(ttheta_seq, -249700 + t_logthetahypprior(ttheta_seq), col = "blue")
legend("topright", legend = c("true value", "unnormalized log half-Cauchy"),
       col = c("red", "blue"), lwd = rep(1, 2))
peaky_gamma <- sapply(ttheta_seq, function(k) integrated_kernel(t_mu = t_mu, t_theta = k, t_pid = t_pid, t_xi = t_xi, D = c(18, 30), A = A, Time = Time, a = a, r = R, e = e))
peaky_avg_gamma <- sapply(1:5,
       function(j) sapply(ttheta_seq, function(k) integrated_kernel(t_mu = emcee_int_thin_flat[j, 1], t_theta = k, t_pid = emcee_int_thin_flat[j, 3], t_xi = emcee_int_thin_flat[j, 4], D = emcee_results_list$simulated_data, A = A, Time = Time, a = a, r = R, e = e)))
instable_index <- sapply(1:5, function(j) which(peaky_gamma > -249700 + t_logthetahypprior(ttheta_seq)))
instable_index <- which(peaky_gamma > -249700 + t_logthetahypprior(ttheta_seq))
#exp(-max(peaky_gamma[-(1:max(instable_index))]) + peaky_gamma[max(instable_index) + 1])
#pcauchy(exp(ttheta_seq[max(instable_index)]), location = 8.7 * 10^-12, scale = 8.7 * 10^-12 * 100, log.p = TRUE) - pcauchy(0, location = 8.7 * 10^-12, scale = 8.7 * 10^-12 * 100, log.p = TRUE)
ptrunc(exp(ttheta_seq[max(instable_index)]), spec = "cauchy", a = 0, location = 8.7 * 10^-12, scale = 8.7 * 10^-12 * 100)

plot(ttheta_seq,
     sapply(ttheta_seq, function(k) {
       lshape <- 2 * t_mu - k
       exp(lshape) * (t_mu - k) - lgamma(exp(lshape)) + logplusvec(lgamma(exp(lshape) + 0:10) + (-exp(t_shape) - 0:10) * (logplus(t_rate, log(R) + log(e) + log(Time))))
     }),
     type = "l",
     xlab = "log theta", ylab = "Laplace(Gamma)", main = "Lazhi parametrized Gamma on theta")
abline(v = t_theta, col = "red")
plot(tshape_seq,
     sapply(tshape_seq, function(k) integrated_kernel(t_shape = k, t_rate = t_rate, t_pid = t_pid, t_xi = t_xi, D = emcee_results_list$simulated_data, A = A, Time = Time, a = a, r = R, e = e)),
     type = "l",
     xlab = "shape parameter alpha", ylab = "marginal posterior kernel", 
     main = "marginal posterior kernel vs. Gamma alpha")
abline(v = t_shape, col = "green")
plot(tshape_seq,
     sapply(tshape_seq, function(k) analytical(k, t_rate, t_pid, t_xi, D = c(0:30, 30), A = A, Time = Time, a = a, r = R, e = e)),
     type = "l",
     xlab = "shape parameter alpha", ylab = "integral", 
     main = "integral vs. Gamma alpha")
abline(v = t_shape, col = "green")
dev.off()
# re-run the hyperprior
pdf("meeting_shared/9Nov2023/final gamma.pdf")
plot(tshape_seq,
     sapply(tshape_seq, function(k) integrated_kernel(t_shape = k, t_rate = t_rate, t_pid = t_pid, t_xi = t_xi, D = c(0:25, 30), A = A, Time = Time, a = a, r = R, e = e)),
     type = "l",
     xlab = "shape parameter alpha", ylab = "marginal posterior kernel", 
     main = "marginal posterior kernel vs. Gamma alpha")
abline(v = t_shape, col = "green")
dev.off()

layout(rbind(c(0, 11, 11, 11, 11),
             c(12, 1, 0, 0, 0),
             c(13, 5, 2, 0, 0),
             c(14, 6, 7, 3, 0),
             c(15, 8, 9, 10, 4),
             c(0, 16, 17, 18, 19)),
       width = c(lcm(2), 1, 1, 1, 1),
       height = c(lcm(2), 1, 1, 1, 1, lcm(2)))
layout.show(19)
layout(rbind(c(0, 11, 11, 11, 11, 11),
             c(12, 0, 1, 20, 20, 0),
             c(13, 0, 5, 2, 0, 0),
             c(14, 0, 6, 8, 3, 0),
             c(15, 0, 7, 9, 10, 4),
             c(0, 0, 16, 17, 18, 19)),
       width = c(lcm(2), lcm(2), 1, 1, 1, 1),
       height = c(lcm(2), 1, 1, 1, 1, lcm(3)))
par(mar = rep(0, 4), las = 1, cex = 1, tck = 0.01)
sapply(1:3, function(i) {
  plot(density(emcee_thin_flat[, i]), main = "", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  abline(v = true_values[i], col = "red")
})
plot(density(emcee_thin_flat[, 4]), main = "", xlab = paste("log", expression(xi)), ylab = "", yaxt = "n")
abline(v = true_values[4], col = "red")
emcee_kde <- kde(emcee_thin_flat[, index_mat[1, ]])
plot(emcee_kde, display = "filled.contour", cont = cont, col.fun = viridis::viridis,
     xlab = "", ylab = paste("log", expression(theta)), xaxt = "n")
plot(emcee_kde, cont = cont, display = "slice", add = TRUE)
emcee_kde <- kde(emcee_thin_flat[, index_mat[2, ]])
plot(emcee_kde, display = "filled.contour", cont = cont, col.fun = viridis::viridis,
     xlab = "", ylab = paste("logit", expression(pi)), xaxt = "n")
plot(emcee_kde, cont = cont, display = "slice", add = TRUE)
emcee_kde <- kde(emcee_thin_flat[, index_mat[3, ]])
plot(emcee_kde, display = "filled.contour", cont = cont, col.fun = viridis::viridis,
     xlab = paste("log", expression(mu)), ylab = paste("log", expression(xi)))
plot(emcee_kde, cont = cont, display = "slice", add = TRUE)
emcee_kde <- kde(emcee_thin_flat[, index_mat[4, ]])
plot(emcee_kde, display = "filled.contour", cont = cont, col.fun = viridis::viridis,
     xlab = "", ylab = "", xaxt = "n", yaxt = "n")
plot(emcee_kde, cont = cont, display = "slice", add = TRUE)
emcee_kde <- kde(emcee_thin_flat[, index_mat[5, ]])
plot(emcee_kde, display = "filled.contour", cont = cont, col.fun = viridis::viridis,
     xlab = paste("log", expression(theta)), ylab = "", yaxt = "n")
plot(emcee_kde, cont = cont, display = "slice", add = TRUE)
emcee_kde <- kde(emcee_thin_flat[, index_mat[6, ]])
plot(emcee_kde, display = "filled.contour", cont = cont, col.fun = viridis::viridis,
     xlab = paste("logit", expression(pi)), ylab = "", yaxt = "n")
plot(emcee_kde, cont = cont, display = "slice", add = TRUE)
#mtext("Density and contour plots of parameters", cex = 1.5, font = 2, outer = TRUE)
pstr = 60
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(0.5, 0.25, "Density and contour plots of parameters", cex = 1.3, font = 2)
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(0.5, 0.5, expression("log("*mu*")"), str = pstr)
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(0.5, 0.5, expression("log("*theta*")"), str = pstr)
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(0.5, 0.5, expression("logit("*pi[d]*")"), str = pstr)
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(0.5, 0.5, expression("log("*xi*")"), str = pstr)
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(0.5, 0.5, expression("log("*mu*")"))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(0.5, 0.5, expression("log("*theta*")"))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(0.5, 0.5,  expression("logit("*pi[d]*")"))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(0.5, 0.5, expression("log("*xi*")"))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(0.5, 0.75, "ordinary model", font = 2, cex = 1.1)
