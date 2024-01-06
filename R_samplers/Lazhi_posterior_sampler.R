library(emg)
library(mvtnorm)
library(MASS)
library(Rcpp)
library(LaplacesDemon)
sourceCpp("logsum.cpp")
source("kernel_functions.R")

Lazhi_simulator <- function(iterations = 1000, t_init = c(t_mu, t_theta, t_pid, t_xi),
                            D, A, Time, a, r, e,
                            alpha0 = 1e-6, beta0 = 1.1, 
                            #n_mh_iter = 100,
                            #emg_lambda = 1e-1,
                            debug_mode = FALSE,
                            transform_parameters = TRUE) {
  n <- length(D) - 1
  init <- c(exp(t_init[1]), exp(t_init[2]), invlogit(t_init[3]), exp(t_init[4]))
  
  #####################
  # internal functions
  ####################
  r0infgamma <- function(n, zeroprob, lmu, ltheta) {
    zeroindexes <- which(log(runif(n)) < log(zeroprob))
    Lambda <- rgamma(n = n, 
                     shape = exp(2 * lmu - ltheta), 
                     rate = exp(lmu - ltheta))
    Lambda[zeroindexes] <- 0
    Lambda
  }
  lQ <- function(lparams, Lambda) {
    nd <- sum(Lambda == 0)
    lmu <- lparams[1]
    ltheta <- lparams[2]
    alpha <- exp(2 * lmu - ltheta)
    lbeta <- lmu - ltheta
    lp <- t_logmuhypprior(lmu) + t_logthetahypprior(ltheta)
    lmid <- (n - nd) * (alpha * lbeta - lgamma(alpha))
    llast <- -exp(lbeta) * sum(Lambda) + (alpha - 1) * sum(log(Lambda[Lambda > 0]))
    lp + lmid + llast - lmu - ltheta # the last two terms are Jacobian
  }
  Q <- function(params, Lambda) {
    if (all(params > 0)) {
      lparams <- log(params)
      lQ(lparams, Lambda)
    } else {-Inf}
  }
  
  ##############
  # initial step (values generated at each step+needed for the next step that are not stored)
  ##############
  Lambda <- r0infgamma(n = n, zeroprob = init[3], lmu = t_init[1], ltheta = t_init[2])
  nd <- rep(NA, iterations + 1)
  nd[1] <- sum(Lambda == 0)
  
  # resulting_sample = stack of vectors (log(mu), log(theta), logit(pid), log(xi)) on the original scale
  resulting_sample <- matrix(nrow = iterations + 1, ncol = 4)
  resulting_sample[1, ] <- t_init
  
  it <- 2
  out_break <- FALSE
  while (it <= (iterations + 1)) {
    #########
    # step 1
    #########
    resulting_sample[it, 3] <- rbeta(n = 1, shape1 = nd + 1, shape2 = n - nd + 1)
    
    #########
    # step 2
    ########
    laxi <- log(a) + resulting_sample[it - 1, 4]
    lomega <- sapply(1:n, function(j) laxi - logplus(laxi, log(r) + log(e) + log(Lambda[j])))
    background_prob <- ifelse(all(lomega != 0), exp(lomega),
                              (laxi * Time)/D[-length(D)])
    Background <- rbinom(n = n, size = D[-length(D)], prob = exp(lomega))
    
    ########
    # step 3
    ########
    alphak <- alpha0 + D[length(D)] + sum(Background)
    betak <- beta0 + A * Time + a * Time
    resulting_sample[it, 4] <- rgamma(n = 1, shape = alphak, rate = betak) |> log()
    
    ########
    # step 4
    ########
    # derivatives
    if (transform_parameters) { # mh on transformed parameter space
      opt_list <- optim(par = t_init[1:2],#resulting_sample[it - 1, 1:2],
                        fn = lQ, method = "BFGS", control = list(fnscale = -1), hessian = TRUE, Lambda = Lambda)
      tau <- if(!out_break) {opt_list$par} else {t_init[1:2]}
      Sigma <- if(!out_break) {-ginv(opt_list$hessian)} else{0.01*diag(nrow = 2, ncol = 2)}
      
      if (!is.symmetric.matrix(Sigma)) Sigma <- LaplacesDemon::as.symmetric.matrix(Sigma)
      if (!is.positive.semidefinite(Sigma)) Sigma <- LaplacesDemon::as.positive.semidefinite(Sigma)
      
      ############################
      # normal-exponantial-modified-Gaussian proposal
      ###########################
      #correlation <- Sigma[1, 2] / sqrt(Sigma[1, 1] * Sigma[2, 2])
      #mu_proposal <- rnorm(n = 1, mean = tau[1], sd = sqrt(Sigma[1, 1]))
      #theta_proposal <- remg(1, mu = tau[2] + correlation * sqrt(Sigma[2, 2]) * (t_mu - tau[1]) / sqrt(Sigma[1, 1]), sigma = sqrt((1 - correlation^2) * Sigma[2, 2]), lambda = emg_lambda)
      #proposal <- c(mu_proposal, theta_proposal)
      #lacceptance_ratio <- lQ(lparams = proposal, Lambda) +
      #  dnorm(x = resulting_sample[it - 1, 1], mean = tau[1], sd = sqrt(Sigma[1, 1]), log = TRUE) +
      #  demg(x = resulting_sample[it - 1, 2], mu = tau[2] + correlation * sqrt(Sigma[2, 2]) * (t_mu - tau[1]) / sqrt(Sigma[1, 1]), sigma = sqrt((1 - correlation^2) * Sigma[2, 2]), lambda = emg_lambda, log = TRUE) - 
      #  lQ(lparams = resulting_sample[it - 1, 1:2], Lambda) -
      #  dnorm(x = mu_proposal, mean = tau[1], sd = sqrt(Sigma[1, 1]), log = TRUE) -
      #  demg(x = theta_proposal, mu = tau[2] + correlation * sqrt(Sigma[2, 2]) * (t_mu - tau[1]) / sqrt(Sigma[1, 1]), sigma = sqrt((1 - correlation^2) * Sigma[2, 2]), lambda = emg_lambda, log = TRUE)
      #lacceptance_prob <- min(0, lacceptance_ratio)
      #accepted <- log(runif(1)) < lacceptance_prob
      
      #if (!is.na(accepted)) {
      #  if (accepted) {
      #    resulting_sample[it, 1:2] <- proposal
      #  } else {
      #    resulting_sample[it, 1:2] <- resulting_sample[it - 1, 1:2]
      #  }
      #} else {
      #  resulting_sample[it, 1:2] <- resulting_sample[it - 1, 1:2]
      #}
      
      #print(list(accept_numlQ = lQ(lparams = proposal, Lambda),
      #           accept_numdnorm = dnorm(x = resulting_sample[it - 1, 1], mean = tau[1], sd = sqrt(Sigma[1, 1]), log = TRUE),
      #           accept_numemg = demg(x = resulting_sample[it - 1, 2], mu = tau[2] + correlation * sqrt(Sigma[2, 2]) * (t_mu - tau[1]) / sqrt(Sigma[1, 1]), sigma = sqrt((1 - correlation^2) * Sigma[2, 2]), lambda = emg_lambda, log = TRUE),
      #           accept_denomlQ = lQ(lparams = resulting_sample[it - 1, 1:2], Lambda),
      #           accept_denomdnorm = dnorm(x = mu_proposal, mean = tau[1], sd = sqrt(Sigma[1, 1]), log = TRUE),
      #           accept_denomemg = demg(x = theta_proposal, mu = tau[2] + correlation * sqrt(Sigma[2, 2]) * (t_mu - tau[1]) / sqrt(Sigma[1, 1]), sigma = sqrt((1 - correlation^2) * Sigma[2, 2]), lambda = emg_lambda, log = TRUE)))
      
      ################################################
      # multiple mh steps
      ################################
      #mh_sample <- matrix(nrow = n_mh_iter + 1, ncol = 2)
      #mh_sample[1, ] <- resulting_sample[it - 1, 1:2]
      #for (o in 2:(n_mh_iter + 1)) {
      #  proposal <- rmvnorm(n = 1, mean = tau, sigma = Sigma)
      #  lacceptance_ratio <- lQ(lparams = mh_sample[o - 1, ], Lambda) +
      #    dmvnorm(x = exp(mh_sample[o - 1, ]), mean = tau, sigma = Sigma, log = TRUE) - 
      #    lQ(lparams = proposal, Lambda) -
      #    dmvnorm(x = exp(proposal), mean = tau, sigma = Sigma, log = TRUE)
      #  lacceptance_prob <- min(0, lacceptance_ratio)
      #  accepted <- log(runif(1)) < lacceptance_prob
      
      #  if (!is.na(accepted)) {
      #    if (accepted) {
      #      mh_sample[o, ] <- proposal
      #    } else {
      #      mh_sample[o, 1:2] <- mh_sample[o - 1, 1:2]
      #    }
      #  } else {
      #    mh_sample[o, 1:2] <- mh_sample[o - 1, 1:2]
      #  }
      #}
      #resulting_sample[it, 1:2] <- mh_sample[sample.int(n_mh_iter, size = 1) + 1, ]
      
      #################################################################
      # classical bivariate normal proposal
      ###########################################
      #lacceptance_ratio <- rep(NA, iterations)
      proposal <- rmvnorm(n = 1, mean = tau, sigma = Sigma)
      lacceptance_ratio <- lQ(lparams = proposal, Lambda) +
        dmvnorm(x = resulting_sample[it - 1, 1:2], mean = tau, sigma = Sigma, log = TRUE) - 
        lQ(lparams = resulting_sample[it - 1, 1:2], Lambda) -
        dmvnorm(x = proposal, mean = tau, sigma = Sigma, log = TRUE)
      lacceptance_prob <- min(0, lacceptance_ratio)
      accepted <- log(runif(1)) < lacceptance_prob
      
      if (!is.na(lacceptance_ratio)) {
        if (accepted) {
          resulting_sample[it, 1:2] <- proposal
        } else {
          resulting_sample[it, 1:2] <- resulting_sample[it - 1, 1:2]
        }
        out_break <- FALSE # has the chain drifted to singularity?
      } else {
        resulting_sample[it, 1:2] <- resulting_sample[it - 1, 1:2]
        warning(paste("chain drifted to singularity (all lambda small) at iteration ", it - 1))
        out_break <- TRUE
      }
      
      ########################################################################
      # exponential proposal (on shape, rate)
      #####################################
      #beta_proposal <- rexp(n = 1, rate = sum(Lambda))
      #alpha_proposal <- 1 - rexp(n = 1, rate = -sum(log(Lambda[Lambda > 0])))
      #proposal <- c(alpha_proposal, beta_proposal)
      #initial_alpha <- exp(2 * resulting_sample[it - 1, 1] - resulting_sample[it - 1, 2])
      #initial_beta <- exp(resulting_sample[it - 1, 1] - resulting_sample[it - 1, 2])
      #lacceptance_ratio <- Q_shaperate(shape = initial_alpha, rate = initial_beta, Lambda) +
      #  dexp(initial_beta, rate = sum(Lambda), log = TRUE) + dexp(1 - initial_alpha, rate = -sum(log(Lambda[Lambda > 0])), log = TRUE) - 
      #  Q_shaperate(shape = alpha_proposal, rate = beta_proposal, Lambda) -
      #  dexp(beta_proposal, rate = sum(Lambda), log = TRUE) - dexp(1 - alpha_proposal, rate = -sum(log(Lambda[Lambda > 0])), log = TRUE)
      #lacceptance_prob <- min(0, lacceptance_ratio)
      #accepted <- log(runif(1)) < lacceptance_prob
      
      #if (!is.na(accepted)) {
      #  if (accepted) {
      #    resulting_sample[it, 1:2] <- log(proposal)
      #  } else {
      #    resulting_sample[it, 1:2] <- resulting_sample[it - 1, 1:2]
      #  }
      #} else {
      #  resulting_sample[it, 1:2] <- resulting_sample[it - 1, 1:2]
      #}
    } else { # mh on original parameter space
      opt_list <- optim(par = init[1:2],#resulting_sample[it - 1, 1:2],
                        fn = Q, method = "SANN", control = list(fnscale = -1), 
                        #lower = rep(.Machine$double.xmin, 2), upper = Inf, 
                        Lambda = Lambda)
      opt_hess <- numDeriv::hessian(func = Q, x)
      tau <- if(!out_break) {opt_list$par} else {init[1:2]}
      Sigma <- if(!out_break) {-ginv(opt_hess)} else{0.01*diag(nrow = 2, ncol = 2)}
      
      if (!is.symmetric.matrix(Sigma)) Sigma <- LaplacesDemon::as.symmetric.matrix(Sigma)
      if (!is.positive.semidefinite(Sigma)) Sigma <- LaplacesDemon::as.positive.semidefinite(Sigma)
      
      proposal <- rmvnorm(n = 1, mean = tau, sigma = Sigma)
      lacceptance_ratio <- Q(params = proposal, Lambda) +
        dmvnorm(x = exp(resulting_sample[it - 1, 1:2]), mean = tau, sigma = Sigma, log = TRUE) - 
        Q(params = exp(resulting_sample[it - 1, 1:2]), Lambda) -
        dmvnorm(x = proposal, mean = tau, sigma = Sigma, log = TRUE)
      lacceptance_prob <- min(0, lacceptance_ratio)
      accepted <- log(runif(1)) < lacceptance_prob
      
      if (!is.na(lacceptance_ratio)) {
        if (accepted) {
          resulting_sample[it, 1:2] <- log(proposal)
        } else {
          resulting_sample[it, 1:2] <- resulting_sample[it - 1, 1:2]
        }
        out_break <- FALSE # has the chain drifted to singularity?
      } else {
        resulting_sample[it, 1:2] <- resulting_sample[it - 1, 1:2]
        warning(paste("chain drifted to singularity (all lambda small) at iteration ", it - 1))
        out_break <- TRUE
      }
    }
    
    
    ########
    # step 5
    ########
    Signal <- D[-length(D)] - Background
    no_signal <- Signal == 0
    len_no_signal <- length(which(Signal == 0))
    
    lalpha_tilt <- 2 * resulting_sample[it, 1] - resulting_sample[it, 2]
    lbeta_tilt <- resulting_sample[it, 1] - resulting_sample[it, 2]
    lalpha_tilti <- sapply(1:length(Signal), function(i) logplus(lalpha_tilt, log(Signal[i])))
    lbeta_tilti <- logplus(lbeta_tilt, log(r) + log(e) + log(Time))
    lpistar_tilt <- log(resulting_sample[it, 3]) - logplus(log(resulting_sample[it, 3]), log(1 - resulting_sample[it, 3]) + exp(lalpha_tilt) * (lbeta_tilt - lbeta_tilti))
    
    Lambda <- rep(NA, n)
    Lambda[no_signal] <- r0infgamma(n = len_no_signal, zeroprob = exp(lpistar_tilt), lmu = lalpha_tilti - lbeta_tilti, ltheta = lalpha_tilti - 2 * lbeta_tilti)
    Lambda[!no_signal] <- rgamma(n = n - len_no_signal, shape = exp(lalpha_tilti), rate = exp(lbeta_tilti))
    
    if(out_break) {
      #it <- max(2, min(max(which(nd <= nd[1])), max(which(resulting_sample[, 3] >= quantile(resulting_sample[, 3], prob = 0.95, na.rm = TRUE))), max(which(resulting_sample[, 2] >= quantile(resulting_sample[, 2], prob = 0.95, na.rm = TRUE)))) - 100)
      it <- max(2, which(nd <= nd[1]))
      Lambda <- r0infgamma(n = n, zeroprob = init[3], lmu = t_init[1], ltheta = t_init[2])
    } else {
      it <- it + 1
    }
    nd[it] <- sum(Lambda == 0)
    
    if (((it - 1) %% 100) == 0) print(paste("iteration", it - 1))
    
    if (debug_mode) {
      print(list(omega = exp(lomega), Lambda = Lambda,
                 Sigma_symmetric = isSymmetric(Sigma),
                 Sigma_positivesemidefinite = is.positive.semidefinite(Sigma),
                 accept_numlQ = lQ(lparams = proposal, Lambda),
                 accept_nummvnorm = dmvnorm(x = resulting_sample[it - 1, 1:2], mean = tau, sigma = Sigma, log = TRUE),
                 accept_denomlQ = lQ(lparams = resulting_sample[it - 1, 1:2], Lambda),
                 accept_denommvnorm = dmvnorm(x = proposal, mean = tau, sigma = Sigma, log = TRUE),
                 initial_state = resulting_sample[it - 1, 1:2], proposal = proposal, acceptance_ratio = exp(lacceptance_ratio), acception = accepted,
                 Signal = Signal, zeroprobability = exp(lpistar_tilt), pid = resulting_sample[it, 3], 
                 lmu = lalpha_tilti - lbeta_tilti, ltheta = lalpha_tilti - 2 * lbeta_tilti, alpha = exp(lalpha_tilt), beta = exp(lbeta_tilt), betai = exp(lbeta_tilti)))
    }
  }
  colnames(resulting_sample) <- c("t_mu", "t_theta", "pid", "t_xi")
  resulting_sample[-1, ]
}
set.seed(49)
L <- Lazhi_simulator(iterations = 100, D = emcee_results_list$simulated_data, A = A, Time = Time, a = a, r = R, e = e, transform_parameters = FALSE)

#t_theta_x <- seq(-25, -5, by = 0.1)
#plot(t_theta_x, exp(sapply(t_theta_x, function(i) lQ(c(t_mu, i), Lambda))), type = "l")
#abline(v = t_theta, col = "red")
#opt_list <- optim(par = t_init[1:2],#resulting_sample[it - 1, 1:2],
#                  fn = lQ, method = "BFGS", control = list(fnscale = -1), hessian = TRUE, Lambda = Lambda)
#tau <- opt_list$par
#Sigma <- -ginv(opt_list$hessian)
#correlation <- Sigma[1, 2] / sqrt(Sigma[1, 1] * Sigma[2, 2])
#lines(t_theta_x, dnorm(t_theta_x, mean = tau[2] + correlation * sqrt(Sigma[2, 2]) * (t_mu - tau[1]) / sqrt(Sigma[1, 1]), sd = sqrt((1 - correlation^2) * Sigma[2, 2]), log = TRUE), col = "green")
#lines(t_theta_x, (demg(t_theta_x, mu = tau[2] + correlation * sqrt(Sigma[2, 2]) * (t_mu - tau[1]) / sqrt(Sigma[1, 1]), sigma = sqrt((1 - correlation^2) * Sigma[2, 2]), lambda = 1e1, log = TRUE)), col = "purple", type = "l")

#plot(seq(-15, -7, by = 0.1), (sapply(seq(-15, -7, by = 0.1), function(i) lQ(c(i, t_theta), Lambda))), type = "l")
#abline(v = t_mu, col = "red")
#lines(t_theta_x, dnorm(t_theta_x, mean = tau[2] + correlation * sqrt(Sigma[2, 2]) * (t_mu - tau[1]) / sqrt(Sigma[1, 1]), sd = sqrt((1 - correlation^2) * Sigma[2, 2]), log = TRUE), col = "green")
#plot(seq(-50, 50, by = 0.1), dnorm(seq(-50, 50, by = 0.1), log = TRUE), type = "l")

#x     <- seq(-21, -5, 0.1) 
#y     <- seq(-15, -7, 0.1)
#f     <- function(x, y) lQ(c(x, y), Lambda)
#z     <- sapply(x, function(j) sapply(y, function(i) f(j, i)))

#create contour plot
#contour(x = y, y = x, z = z, xlab = "t_theta", ylab = "t_mu")
#abline(h = t_mu, col = "red")
#abline(v = t_theta, col = "red")

#LaplacesDemon(Model = function(parm, Data) lQ(lparams = parm, Lambda = Data), Data = Lambda, Initial.Values = t_init[1:2], Covar = Sigma, Iterations = 11, Thinning = 2, Algorithm = "SGLD")
