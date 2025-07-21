
library(mvtnorm)
library(numDeriv)
library(matrixcalc)
library(tictoc)
library(scales)
library(ggplot2)
library(GGally)
library(patchwork)
#library(reticulate)


source("functions.R")
source("data.simulate.R")
source("gibbs-new.R")
source("data.xray.R")
source("data.optical.R")
#source("design.response.R")

Data <- simulate.data(pop.mean =10, pop.sd=5, pop.zero.prob=0.4, n.segments=200)

david <- FALSE
if (david) {
  xData <- readRDS("xData.RDS")
  oData <- readRDS("oData.RDS")
} else {
  xData <- xray.data(count_file = "/Users/jl1023/NGC2516/27/repro/xYroi.txt",
                     area_file = "/Users/jl1023/NGC2516/27/repro/xAroi.txt",
                     exposure_file = "/Users/jl1023/NGC2516/27/repro/xEroi.txt",
                     rmat_file = "/Users/jl1023/NGC2516/27/repro/xRmat.txt",
                     bkgd_file = "/Users/jl1023/NGC2516/27/cochran/cbind_bkgd.txt")
  
  oData <- optical.data(count_file = "/Users/jl1023/NGC2516/27/cochran/optYroi.txt",
                        area_file = "/Users/jl1023/NGC2516/27/cochran/optAroi.txt",
                        exposure_file = "/Users/jl1023/NGC2516/27/cochran/optEroi.txt",
                        rmat_npz = "/Users/jl1023/NGC2516/27/cochran/optRmat.txt.npz",
                        bkgd_file = "/Users/jl1023/NGC2516/27/cochran/cbind_bkgd.txt")
  
}

# sampler becomes unstable if population is too concentrated near zero
set.seed(43)
tic("Xray data")
xMCMC.sample <- gibbs(seg=xData$Segments, 
                     sources=xData$Sources, 
                     bkgd=xData$Background)
toc(quiet = FALSE, log = TRUE)
(log.txt <- tic.log(format = TRUE))
(log.lst <- tic.log(format = TRUE))
tic.clearlog()

set.seed(42)
tic("optical data")
oMCMC.sample <- gibbs(seg=oData$Segments, 
                      sources=oData$Sources, 
                      bkgd=oData$Background)
toc(quiet = FALSE, log = TRUE)
(log.txt <- tic.log(format = TRUE))
(log.lst <- tic.log(format = TRUE))
tic.clearlog()

pdf(file = "traceXrayPosterior.pdf")
plot(xMCMC.sample$pop.shape[-(1:100)], type = "l", ylab = "population shape", main = "X-ray posterior traceplot")
plot(xMCMC.sample$bkgd.rate[-(1:100)], type = "l", ylab = "background rate", main = "X-ray posterior traceplot")
plot(xMCMC.sample$pop.rate[-(1:100)], type = "l", ylab = "population rate", main = "X-ray posterior traceplot")
plot(xMCMC.sample$pop.zero.prob[-(1:100)], type = "l", ylab = "population 0 prob", main = "X-ray posterior traceplot")
#abline(h = 1-86/131, col = "red")
plot(xMCMC.sample$pop.mean[-(1:100)], type = "l", ylab = "population mean", main = "X-ray posterior traceplot")
plot(xMCMC.sample$pop.sd[-(1:100)], type = "l", ylab = "population sd", main = "X-ray posterior traceplot")
plot(xMCMC.sample$pop.mean[-(1:100)], xMCMC.sample$pop.sd[-(1:100)], type = "l", xlab = "population mean", ylab = "population sd")
for (i in 1:xData$Sources$n) {
  plot(xMCMC.sample$src.rates[-(1:100), i], type = "l", ylab = paste("source intensity", i), main = "X-ray posterior traceplot")
}
x0sources <- colSums(xMCMC.sample$zero) / 5e4
hist(x0sources, breaks = seq(0, max(x0sources)+0.01, by = 0.01), 
     xlab = "proportion of iterations each lambda is 0", 
     main = "Proportion of dark iterations",
     sub = "for each X-ray detection")
dev.off()

pdf(file = "traceOpticalPosterior.pdf")
plot(oMCMC.sample$pop.shape[-(1:100)], type = "l", ylab = "population shape", main = "optical posterior traceplot")
plot(oMCMC.sample$bkgd.rate[-(1:100)], type = "l", ylab = "background rate", main = "optical posterior traceplot")
plot(oMCMC.sample$pop.rate[-(1:100)], type = "l", ylab = "population rate", main = "optical posterior traceplot")
plot(oMCMC.sample$pop.zero.prob[-(1:100)], type = "l", ylab = "population 0 prob", main = "optical posterior traceplot")
abline(h = quantile(oMCMC.sample$pop.zero.prob[-(1:100)], probs = 0.025), col = "red")
abline(h = quantile(oMCMC.sample$pop.zero.prob[-(1:100)], probs = 0.975), col = "red")
plot(oMCMC.sample$pop.mean[-(1:100)], type = "l", ylab = "population mean", main = "optical posterior traceplot")
plot(oMCMC.sample$pop.sd[-(1:100)], type = "l", ylab = "population sd", main = "optical posterior traceplot")
plot(oMCMC.sample$pop.mean[-(1:100)], oMCMC.sample$pop.sd[-(1:100)], type = "l", xlab = "population mean", ylab = "population sd")
for (i in 1:10) {
  plot(oMCMC.sample$src.rates[-(1:100), i], type = "l", ylab = paste("source intensity", i), main = "optical  posterior traceplot")
}

o0sources <- colSums(oMCMC.sample$zero) / (nrow(oMCMC.sample$src.rates) - 1)
hist(o0sources, breaks = seq(0, 1, by = 0.05), 
     xlab = "proportion of iterations each lambda is 0", 
     main = "Proportion of dark iterations",
     sub = "for each optical source")

plot(oMCMC.sample$src.rates[-(1:100), 1190], type = "l", ylab = "source intensity 1190 and 1202", 
     main = "optical  posterior traceplot", sub = "non-zero optical source intensities")
lines(oMCMC.sample$src.rates[-(1:100), 1202], col = "red")

plot(rowSums(oMCMC.sample$src.rates)[-(1:100)], type = "l", main = "Aggregated traceplot of all sources", ylab = "summed intensities")
plot(rowSums(oMCMC.sample$src.rates[, o0sources < 0.9])[-(1:100)], type = "l", main = "Aggregated traceplot of bright sources", ylab = "summed intensities")
plot(rowSums(oMCMC.sample$src.rates[, o0sources > 0.9])[-(1:100)], type = "l", main = "Aggregated traceplot of dark sources", ylab = "summed intensities")
plot(rowSums(oMCMC.sample$src.rates[, o0sources > 0.9])[-(1:100)] / rowSums(oMCMC.sample$src.rates)[-(1:100)], type = "l", main = "Aggregated traceplot of dark sources", ylab = "intensity proportion")
o0wrapper <- function(thre, small = TRUE) {
  if (small) {
    sapply(thre, FUN = function(i) {mean(rowSums(oMCMC.sample$src.rates[, o0sources <= i])[-(1:100)])})
  } else {
    sapply(thre, FUN = function(i) {mean(rowSums(oMCMC.sample$src.rates[, o0sources > i])[-(1:100)])})
  }
}
eval.seq <- seq(0, 1, by = 0.01)
plot(eval.seq, o0wrapper(eval.seq, small = TRUE) / o0wrapper(1, small = TRUE), type = "l",
     xlab = "proportion of non-0 iterations", ylab = "intensities proportion", main = "Proportion of mean aggregated intensities", sub = "for filtering segments")

singleSource <- oData$Segments$cnt.obs[(oData$Segments$src.index %in% which(o0sources < 0.9))[oData$Segments$n.overlap == 1]]
single09Exposure <- oData$Segments$eff.area[(oData$Segments$src.index %in% which(o0sources < 0.9))[oData$Segments$n.overlap == 1]]
single09Lumi <- oData$Segments$r[, 1][(oData$Segments$src.index %in% which(o0sources < 0.9))[oData$Segments$n.overlap == 1]]
single09Rate <- singleSource / (single09Exposure * single09Lumi)
single09Area <- oData$Segments$area[(oData$Segments$src.index %in% which(o0sources < 0.9))[oData$Segments$n.overlap == 1]]

plot(oMCMC.sample$pop.sd[-(1:100)], oMCMC.sample$pop.zero.prob[-(1:100)], type = "l", 
     xlab = "population sd", ylab = "population 0 prob")
hist(oData$Segments$cnt.obs[oData$Segments$n.overlap == 1], 
     main = "Photons per single-source segment", xlab = "number of photons")
hist(singleSource, 
     add = TRUE, col = "red", breaks = seq(0, 300, by = 20))
legend("topright", col = "red", pch = 19, legend = "bright intensities")
bin.rates <- table(cut(singleSource, breaks = seq(0, 300, by = 20))) / table(cut(oData$Segments$cnt.obs[oData$Segments$n.overlap == 1], breaks = seq(0, 300, by = 20)))
plot(bin.rates, ylab = "proportion of segments", xlab = "number of photons", main = "Bright intesities in each binning", sub = "among single-source segments")

singleIndex <- oData$Segments$src.index[oData$Segments$n.overlap == 1, 1]
single0prop <- (1-o0sources)[singleIndex]
single09Index <- oData$Segments$src.index[(oData$Segments$src.index %in% which(o0sources < 0.9))[oData$Segments$n.overlap == 1]]
single09prop <- (1-o0sources)[single09Index[single09Index != 0]]
singleCounts <- oData$Segments$cnt.obs[oData$Segments$n.overlap == 1]
singleExposure <- oData$Segments$eff.area[oData$Segments$n.overlap == 1]
singleArea <- oData$Segments$area[oData$Segments$n.overlap == 1]
singleLumi <- oData$Segments$r[oData$Segments$n.overlap == 1, 1]
singleRate <- singleCounts / (singleExposure * singleLumi)
################################
# scaled regression diagnostic
###############################
plot(singleRate, single0prop, pch = 20,
     xlab = "Counts / (area*exposure)", ylab = "proportion of non-0 iterations", main = "Iteration proportions vs. segment covariates", sub = "for single-source segments", col = alpha("grey10", alpha = 0.2))
#points(single09Rate, single09prop, col = alpha("red", 0.7))
glm0prop <- glm(single0prop ~ singleRate, family = quasibinomial(link = "logit"))
summary(glm0prop)
#lines(sort(singleCounts), sort(fitted(glm0prop)), type = "l", col = "blue")
#legend("right", col = c("red", "blue"), pch = 1, legend = c("bright intensities", "fitted values"))

singledf <- as.data.frame(cbind(singleCounts, singleExposure, singleArea, singleLumi), row.names = c("counts", "exposure", "area", "photon proportion"))
ggpairs(singledf, mapping = aes(alpha = 0.01))

############################################
# correctly confounded regression diagnostic
############################################
fit0prop <- glm(single0prop ~ singleCounts + singleExposure + singleLumi, family = quasibinomial(link = "logit"))
summary(fit0prop)

#plot(singleExposure, single0prop, pch = 20,xlab = "exposure", ylab = "proportion of non-0 iterations", main = "Iteration proportions vs. segment exposure", sub = "for single-source segments", col = alpha("grey10", alpha = 0.2))
#points(single09Exposure, single09prop, col = alpha("red", 0.5))
#points(singleExposure, fitted(fit0prop), col = alpha("blue", 0.5))
#legend("topright", inset=c(-1.1,0), col = c("red", "blue"), pch = 1, legend = c("bright intensities", "fitted values"), border = NULL)

##################################
# background regression diagnostic
##################################
mod0prop <- glm(single0prop ~ singleArea + singleLumi, family = quasibinomial(link = "logit"))
summary(mod0prop)


df0prop <- as.data.frame(cbind(single0prop, singleCounts, singleExposure, singleLumi, singleArea, fitted(fit0prop), fitted(mod0prop)))
colnames(df0prop) <- c("non0prop", "counts", "exposure", "prop", "area", "signalFit", "bkgdFit")
p1 <- ggplot(data = df0prop, aes(x = counts, y = non0prop)) +
  geom_point(aes(alpha = 0.1, colour = "true")) +
  geom_point(aes(x = counts, y = signalFit, alpha = 0.1, colour = "fitted")) +
  ylab("Proportion of non 0 iterations")

p2 <- ggplot(data = df0prop, aes(x = exposure, y = non0prop)) +
  geom_point(aes(alpha = 0.1, colour = "true")) +
  geom_point(aes(x = exposure, y = signalFit, alpha = 0.1, colour = "fitted")) +
  ylab("Proportion of non 0 iterations")

p3 <- ggplot(data = df0prop, aes(x = prop, y = non0prop)) +
  geom_point(aes(alpha = 0.1, colour = "true")) +
  geom_point(aes(x = prop, y = signalFit, alpha = 0.1, colour = "fitted")) +
  ylab("Proportion of non 0 iterations") +
  xlab("Photon proportion from source")

p4 <- ggplot(data = df0prop, aes(x = area, y = non0prop)) +
  geom_point(aes(alpha = 0.1, colour = "true")) +
  geom_point(aes(x = area, y = bkgdFit, alpha = 0.1, colour = "fitted")) +
  ylab("Proportion of non 0 iterations")

p1 + p2 + p3 + p4 + 
  plot_layout(guides = 'collect', axes = 'collect') +
  plot_annotation("0 Iteration proportion in sampler vs. photon observations feed in")
dev.off()

saveRDS(list(xMCMC.sample, oMCMC.sample), "MCMCsamples.RDS")
#MCMCsamples <- readRDS("MCMCsamples.RDS"); xMCMC.sample <- MCMCsamples[[1]]; oMCMC.sample <- MCMCsamples[[2]]

for (i in which(o0sources < 0.05)) {
  plot(oMCMC.sample$src.rates[-(1:100), i], type = "l", ylab = paste("source intensity", i), 
       main = "optical  posterior traceplot", sub = "non-zero optical source intensities")
  plot(oMCMC.sample$src.rates[-(1:100), i], oMCMC.sample$pop.mean[-(1:100)], type = "l", col = alpha("grey10", alpha = 0.1),
       xlab = paste0("soure intensity ", i), ylab = "population mean", main = paste0("Source ", i, " vs population mean"))
  plot(oMCMC.sample$src.rates[-(1:100), i], oMCMC.sample$pop.sd[-(1:100)], type = "l", col = alpha("grey10", alpha = 0.1),
       xlab = paste0("soure intensity ", i), ylab = "population sd", main = paste0("Source ", i, " vs population sd"))
}
fitcount <- glm(singleCounts ~ singleExposure + singleLumi * singleArea, family = quasipoisson())
summary(fitcount)
