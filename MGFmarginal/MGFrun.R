library(tictoc)
library(reticulate)
library(Rcpp)

source("data.xray.R")
source("data.optical.R")
source("design.response.R")
source("mgf.R")
sourceCpp("logsum.cpp")

xDesignRes <- xray.design.response(design_file = "/Users/jl1023/NGC2516/27/repro/xDesign.txt",
                                   response_file = "/Users/jl1023/NGC2516/27/repro/xYroi.txt")
oDesignRes <- optical.design.response(design_file = "/Users/jl1023/NGC2516/27/cochran/optDesign.txt.npz",
                                      response_file = "/Users/jl1023/NGC2516/27/cochran/optYroi.txt")

# Computation and verification on solated sources
#+++++++++++++++++++++
#+ NegBin verification
#+++++++++++++++++++++
sourceShape <- (1.2e-6)^2 / 8.7e-12
sourceRates <- 1.2e-6/8.7e-12
sourceProbs <- sourceRates / (sourceRates + xData$Segments$Time * 0.9 * xSData$Isolate$eff.rea)
bkgdShape <- 1e-06 + xData$Background$cnt
bkgdRates <- xData$Segments$Time * xData$Background$area * (1+1e-12)
bkgdProbs <- bkgdRates / (xData$Segments$Time * xSData$Isolate$area + bkgdRates)

tic()
lmgfMarg(xSData, alpha = sourceShape, beta = sourceRates, Time = xData$Segments$Time, bkgdalpha = bkgdShape, bkgdbeta = bkgdRates)
toc()
tic.clear()
tic()
# convoluted negative-binomial true value
sapply(1:length(xSData$Isolate$cnt.obs), 
       FUN = function(i) logplusvec(dnbinom(0:xSData$Isolate$cnt.obs[i], size = sourceShape, prob = sourceProbs[i], log = TRUE) + dnbinom(xSData$Isolate$cnt.obs[i]:0, size = bkgdShape, prob = bkgdProbs[i], log = TRUE))) |> sum()
toc()

# examplar symbolic derivative
eval({xx <- -3; D(expression(xx^3 + xx^2), "xx")})
