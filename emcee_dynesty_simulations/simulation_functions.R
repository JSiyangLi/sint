library(matrixStats)
library(LaplacesDemon) # for logistic function
library(VGAM) # truncated Pareto distribution
library(Rcpp)
sourceCpp("logsum.cpp")

#######################
# simulation 3.2 mh
###########
# static simulation settings
Time <- 5e4
n <- 10
break_point <- 5

########################
# gamma lambdas
Lambdastar_simulator <- function(n, pid, mu_star = 15, theta_star) {
  ########
  # step 2
  ########
  zeroindexes <- which(log(runif(n)) < log(pid))
  Lambda_stars <- rgamma(n, 
                         shape = mu_star^2 / theta_star, 
                         rate = mu_star / theta_star)
  Lambda_stars[zeroindexes] <- 0
  
  Lambda_stars
}

simulator <- function(n, Lambda_stars, xi) {
  ########
  # step 1
  ########
  X <- rpois(1, A * xi * Time)
  
  ########
  # step 3
  ########
  Background <- rpois(n, a * xi * Time)
  Signal <- rpois(n, Lambda_stars)
  Y <- Background + Signal
  
  c(Y, X)
}

##############################
# dynamic simulation settings
##############################
# gamma
######################################
# initial and true values
A <- 2.5e7
xi <- 2e-7
pid = 0.5
theta_star = 50
xi_star = 15
mu_star = 15
a = unlist(xi_star * 100) # 100 from true values of Time and A
R = 1 # modify for estimating Lambda, mu, theta
e = 1 # modify for estimating Lambda, mu, theta
mu = mu_star / (R * e * Time) # true value
theta = theta_star / (R * e * Time)^2 # true value
shape = mu^2 / theta
rate = mu / theta

# transform initial values
t_mu = log(mu)
t_theta = log(theta)
t_pid = LaplacesDemon::logit(pid)
t_xi = log(xi)
t_shape = log(shape)
t_rate = log(rate)

##################################
# NEW! dynamic simulation settings for overlapping and nonhomogeneous background
#############################
# initial and true values
dimI <- 8 # over-writes the previous n
n <- 15
AK = c(22029408, 14093856, 26448800)
X = c(219962, 146332, 285300)
K <- length(AK) # number of background grids
xi_mle = X / (AK * Time)
t_xi_vec_K <- log(xi_mle)
aS = c(700, 200, 1000, 500, 100, 500, 200, 200, 700, 600, 500, 400, 1100, 1500, 1500)
m <- rbind(c(rep(1, 13), 0, 0),
           c(rep(0, 13), 1, 0),
           c(rep(0, 13), 0, 1))
RI <- rbind(c(0.8, 0.05, 0, 0.1, 0.05, rep(0, 10)),
            c(rep(0, 3), 0.2, 0.1, 0.2, 0.2, 0.3, rep(0, 7)),
            c(rep(0, 7), 0.3, 0.3, 0.4, rep(0, 5)),
            c(0, 0.1, 0.75, 0, 0.05, 0, 0.1, rep(0, 8)),
            c(rep(0, 9), 0.5, 0.25, 0.25, rep(0, 3)),
            c(rep(0, 11), 0.4, 0.6, rep(0, 2)),
            c(rep(0, 13), 1, 0),
            c(rep(0, 14), 1)) |> t()
eS <- rep(c(0.9, 1, 1.1), 5)
# no xi_star_over used (refer to xi_mle), mu_star_over or theta_star_over (not needed)

# Lambdastar generator for overlapping by collapsing columns of re
Lambdastar_over_simulator <- function(n, I, eS, RI, pid, mu_star, theta_star) {
  Lambda <- Lambdastar_simulator(I, pid, mu_star, theta_star) / (R * e)
  c(rowSums(eS %*% t(rep(1, n))) / n * (RI %*% Lambda)) # otherwise, rowSums(eS %*% t(rep(1, n)) * RI) forms a vector
}
overlap_simulator <- function(n, Lambda_stars, xi_vec, m_mat, X) {
  ########
  # step 1
  ########
  # neglected since X's are observed in Lazhi's notes.
  
  ########
  # step 3
  ########
  Back_rate <- (rep(1, 3) %*% t(aS)) * m_mat * (xi_vec %*% t(rep(1, 15)))
  Background <- rpois(n, Back_rate[Back_rate > 0] * Time)
  Signal <- rpois(n, Lambda_stars)
  Y <- Background + Signal
  
  Y
}
L <- Lambdastar_over_simulator(n, dimI, eS, RI, pid, mu_star, theta_star)
overlap_simulator(n = 15, L, xi_mle, m_mat = m, X = X)

