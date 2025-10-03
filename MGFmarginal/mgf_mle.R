library(tictoc)
library(Rcpp)
library(LaplacesDemon)

source("data.xray.R")
source("design.response.R")
source("mgf.R")
sourceCpp("logsum.cpp")

#####################
# likelihood wrappers
####################
# log likelihood defined with log parameter space
loglik_lme <- function(lparams, SData, bkgdalpha, bkgdbeta, Time) {
  lalpha <- lparams[1]
  lbeta <- lparams[2]
  if (
    any(c(lalpha, lbeta) <= -.Machine$double.xmax) ||
    any(c(lalpha, lbeta) >= .Machine$double.xmax)
  ) {
    return(as.numeric(-1e100))
  } else {
    l <- lmgfMarg(
      SData,
      lalpha = lalpha,
      lbeta = lbeta,
      Time = Time,
      bkgdalpha = bkgdalpha,
      bkgdbeta = bkgdbeta
    )
    return(l)
  }
}
loglik_lme_meanvar <- function(lparams, SData, bkgdalpha, bkgdbeta, Time) {
  lmu <- lparams[1]
  ltheta <- lparams[2]
  lalpha <- 2 * lmu - ltheta
  lbeta <- lmu - ltheta
  if (
    any(c(lalpha, lbeta) <= -.Machine$double.xmax) ||
    any(c(lalpha, lbeta) >= .Machine$double.xmax)
  ) {
    return(as.numeric(-1e100))
  } else {
    l <- lmgfMarg(
      SData,
      lalpha = lalpha,
      lbeta = lbeta,
      Time = Time,
      bkgdalpha = bkgdalpha,
      bkgdbeta = bkgdbeta
    )
    return(l)
  }
}


############
# X-ray data
############

# xray data, in Gibbs format
xData <- xray.data(
  count_file = "/Users/jl1023/NGC2516/27/repro/xYroi.txt",
  area_file = "/Users/jl1023/NGC2516/27/repro/xAroi.txt",
  exposure_file = "/Users/jl1023/NGC2516/27/repro/xEroi.txt",
  rmat_file = "/Users/jl1023/NGC2516/27/repro/xRmat.txt",
  bkgd_file = "/Users/jl1023/NGC2516/27/cochran/cbind_bkgd.txt"
)

xSData <- matrixSplit(xData)
all(xSData$Isolate$r[, 1] == 0.9) # the density should peak around 0.9


# maximum likelihood
(xlgammamle <- optim(par = rep(0, 2), fn = loglik_lme, 
                     SData = xSData, bkgdalpha = bkgdShape, bkgdbeta = bkgdRates, Time = xData$Segments$Time,
      control = list(fnscale = -1), hessian = TRUE))
(xlgamma_mean <- xlgammamle$par[1] - xlgammamle$par[2]) # mean
(xlgamma_var <- xlgammamle$par[1] - 2 * xlgammamle$par[2]) # variance
(gamma_sd <- sqrt(exp(xlgammamle$par)[1] / exp(xlgammamle$par)[2]^2)) # sd

exp(xlgammamle$par) # alpha and beta
exp(xlgamma_mean)


#################
# optical catalog
#################
# optical data, in Gibbs format
oData <- optical.data(count_file = "/Users/jl1023/NGC2516/27/cochran/optYroi.txt",
                      area_file = "/Users/jl1023/NGC2516/27/cochran/optAroi.txt",
                      exposure_file = "/Users/jl1023/NGC2516/27/cochran/optEroi.txt",
                      rmat_npz = "/Users/jl1023/NGC2516/27/cochran/optRmat.txt.npz",
                      bkgd_file = "/Users/jl1023/NGC2516/27/cochran/cbind_bkgd.txt")
oSData <- matrixSplit(oData)



# maximum likelihood
(olgammamle <- optim(par = rep(0, 2), fn = loglik_lme, 
                     SData = oSData, bkgdalpha = bkgdShape, bkgdbeta = bkgdRates, Time = oData$Segments$Time,
                     control = list(fnscale = -1), hessian = TRUE))
(olgamma_mean <- olgammamle$par[1] - olgammamle$par[2]) # mean
(olgamma_var <- olgammamle$par[1] - 2 * olgammamle$par[2]) # variance
(gamma_sd <- sqrt(exp(olgammamle$par)[1] / exp(olgammamle$par)[2]^2)) # sd

exp(olgammamle$par) # alpha and beta
exp(olgamma_mean)



