library(tictoc)
library(reticulate)
library(Rcpp)

source("data.xray.R")
source("data.optical.R")
source("design.response.R")
sourceCpp("logsum.cpp")

xDesignRes <- xray.design.response(design_file = "/Users/jl1023/NGC2516/27/repro/xDesign.txt",
                                   response_file = "/Users/jl1023/NGC2516/27/repro/xYroi.txt")
oDesignRes <- optical.design.response(design_file = "/Users/jl1023/NGC2516/27/cochran/optDesign.txt.npz",
                                      response_file = "/Users/jl1023/NGC2516/27/cochran/optYroi.txt")

####################
# plotting the background population
####################
MCMCsamples <- readRDS("MCMCsamples.RDS"); xMCMC.sample <- MCMCsamples[[1]]; oMCMC.sample <- MCMCsamples[[2]]

bshape <- 1e-06 + oData$Background$cnt
brate <- oData$Segments$Time * oData$Background$area * (1+1e-12)
bpopmean <- bshape / brate
bvec <- seq(bpopmean - 1e-9, bpopmean + 1e-9, length.out = 1001)

bxshape <- 1e-06 + 1725843
bxrate <- xData$Segments$Time * 160875227.831018 * (1+1e-12)
bxpopmean <- bxshape / bxrate
bxvec <- seq(bxpopmean - 1e-9, bxpopmean + 1e-9, length.out = 1001)

paro <- par()
# Add extra space to right of plot area; change clipping to figure
par(mar=c(10, 4, 4, 1), xpd=TRUE)
plot(density(oMCMC.sample$bkgd.rate[-(1:100)]), xlim = c(min(c(bvec, bxvec)), max(c(bvec, bxvec))),
     xlab = expression(xi), main = "Gamma background posteriors")
lines(density(xMCMC.sample$bkgd.rate[-(1:100)]), #xlim = c(bkgdpopmeanx - 4e-8, bkgdpopmean + 2e-8),
      xlab = expression(xi), main = "Gamma background posteriors", col = "red")
lines(bxvec, dgamma(bxvec, shape = bxshape, rate = bxrate), col = "purple")
lines(bvec, dgamma(bvec, shape = bshape, rate = brate), col = "orange")
legend("topright", col = c("red", "black", "purple", "orange"), lty = 1, inset=c(0,1.7),
       legend = c("X and Y, Xray", "X and Y, optical", "X only, Xray", "X only, Xray&optical"))


par(paro)