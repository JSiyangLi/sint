library(Rcpp)
sourceCpp("logsum.cpp")

source("data.xray.R")
source("data.optical.R")
source("design.response.R")
source("mgf.R")

#=================
# helper functions
#=================
loglik_lme.isolated <- function(lparams, SData) {
  lalpha <- lparams[1]
  lbeta <- lparams[2]
  if (
    any(c(lalpha, lbeta) <= -.Machine$double.xmax) ||
    any(c(lalpha, lbeta) >= .Machine$double.xmax)
  ) {
    return(as.numeric(-1e100))
  } else {
    bkgdShape <- 1e-06 + xData$Background$cnt
    bkgdRates <- xData$Segments$Time * xData$Background$area * (1 + 1e-12)
    l <- lmgfMarg.isolated(
      SData,
      lalpha = lalpha,
      lbeta = lbeta,
      Time = oData$Segments$Time,
      bkgdalpha = bkgdShape,
      bkgdbeta = bkgdRates
    )
    return(l)
  }
}

writeAmat <- function(Glist) {
  # Get all unique segments across all sources
  all_segments <- unique(unlist(lapply(Glist, function(x) x[, "segment"])))
  segments <- sort(all_segments)
  
  # Get all source names
  sources <- names(Glist)
  
  # Create a mapping of segment to e value (take the first occurrence)
  e_values <- numeric()
  for (source in sources) {
    source_data <- Glist[[source]]
    for (i in 1:nrow(source_data)) {
      segment <- as.character(source_data[i, "segment"])
      if (!segment %in% names(e_values)) {
        e_values[segment] <- source_data[i, "e"]
      }
    }
  }
  
  # Create empty matrix with appropriate dimensions
  Amat <- matrix(0, nrow = length(segments), ncol = length(sources),
                 dimnames = list(segments, sources))
  
  # Fill the matrix with r * e values
  for (source in sources) {
    source_data <- Glist[[source]]
    for (i in 1:nrow(source_data)) {
      segment <- as.character(source_data[i, "segment"])
      r_value <- source_data[i, "r"]
      e_value <- e_values[segment]
      Amat[segment, source] <- r_value * e_value
    }
  }
  
  return(Amat)
}

bkgdAmat <- function(Glist, iterations) {
  # Get all unique segments across all sources and sort them
  all_segments <- sort(unique(unlist(lapply(Glist, function(x) x[, "segment"]))))
  
  # Extract a values for each segment (take the first occurrence from any source)
  a_values <- numeric()
  for (segment in all_segments) {
    # Find the first occurrence of this segment in any source
    for (source in names(Glist)) {
      source_data <- Glist[[source]]
      segment_idx <- which(source_data[, "segment"] == segment)
      if (length(segment_idx) > 0) {
        a_values[as.character(segment)] <- source_data[segment_idx[1], "a"]
        break
      }
    }
  }
  
  # Create the a vector in the correct segment order
  a_vector <- a_values[as.character(all_segments)]
  
  # Create matrix with identical rows repeated 'iterations' times
  Amat <- matrix(rep(a_vector, each = iterations), 
                 nrow = iterations, 
                 ncol = length(a_vector),
                 byrow = FALSE,
                 dimnames = list(NULL, as.character(all_segments)))
  
  return(Amat)
}

extractYvec <- function(Glist) {
  # Get all unique segments across all sources and sort them
  all_segments <- sort(unique(unlist(lapply(Glist, function(x) x[, "segment"]))))
  
  # Extract Y values for each segment (take the first occurrence from any source)
  Y_values <- numeric()
  for (segment in all_segments) {
    # Find the first occurrence of this segment in any source
    for (source in names(Glist)) {
      source_data <- Glist[[source]]
      segment_idx <- which(source_data[, "segment"] == segment)
      if (length(segment_idx) > 0) {
        Y_values[as.character(segment)] <- source_data[segment_idx[1], "Y"]
        break
      }
    }
  }
  
  # Create the Y vector in the correct segment order
  Y_vector <- Y_values[as.character(all_segments)]
  names(Y_vector) <- all_segments
  
  return(Y_vector)
}

#=======================
# observation processing
#=======================
oData <- optical.data(count_file = "/Users/jl1023/NGC2516/27/cochran/optYroi.txt",
                      area_file = "/Users/jl1023/NGC2516/27/cochran/optAroi.txt",
                      exposure_file = "/Users/jl1023/NGC2516/27/cochran/optEroi.txt",
                      rmat_npz = "/Users/jl1023/NGC2516/27/cochran/optRmat.txt.npz",
                      bkgd_file = "/Users/jl1023/NGC2516/27/cochran/cbind_bkgd.txt")

# overlapping data
oSData <- matrixSplit(oData)
oGroup <- indexSplit(oSData$Multiple$src.index)
oGData <- yearSplit(oGroup, oSData)
#testGlist <- summarise_GData(oGData[[6]])
bkgdShape <- 1e-06 + oData$Background$cnt
bkgdRates <- oData$Segments$Time * oData$Background$area * (1+1e-12)

# binary overlappings
length_list <- lapply(oGroup$row_indices, length)
binary_indices <- which(unlist(length_list) == 3)

# MLE
# maximum likelihood
(olgammamle <- optim(par = rep(0, 2), fn = loglik_lme.isolated, SData = oSData,
                    control = list(fnscale = -1), hessian = TRUE))
(olgamma_mean <- olgammamle$par[1] - olgammamle$par[2]) # mean
(olgamma_var <- olgammamle$par[1] - 2 * olgammamle$par[2]) # variance
(ogamma_sd <- sqrt(exp(olgammamle$par)[1] / exp(olgammamle$par)[2]^2)) # sd

beta  <- exp(olgammamle$par[2])
beta2 <- bkgdRates
alpha <- exp(olgammamle$par[1])
alpha2 <- bkgdShape

set.seed(42)
iterations <- 1e6
Nseg <- 3
Nsrc <- 2
inCI <- boundarySim <- boundaryCI <- rep(NA, length(binary_indices))
# iterating over binary overlappings
for (i in seq_along(binary_indices)) {
  cat("iteration: ", i, "\n")
  currentG <- summarise_GData(oGData[[binary_indices[i]]])
  
  #--------------------
  # extracting elements
  #--------------------
  # Extract the two elements from current Glist
  elem1 <- currentG[[1]]
  elem2 <- currentG[[2]]
  
  # Identify overlapping segment (common segment between both elements)
  seg1 <- elem1[, "segment"]
  seg2 <- elem2[, "segment"]
  overlapping_segment <- intersect(seg1, seg2)
  
  if (length(overlapping_segment) != 1) {
    stop("Expected exactly one overlapping segment between the two Glist elements")
  }
  
  # Find indices for overlapping and non-overlapping segments
  overlap_idx1 <- which(seg1 == overlapping_segment)
  nonoverlap_idx1 <- which(seg1 != overlapping_segment)
  
  overlap_idx2 <- which(seg2 == overlapping_segment)
  nonoverlap_idx2 <- which(seg2 != overlapping_segment)
  
  # Extract parameters according to the specified pattern
  yu = elem1[nonoverlap_idx1, "Y"]      # Y from non-overlapping segment of first element
  yv = elem2[nonoverlap_idx2, "Y"]      # Y from non-overlapping segment of second element  
  yw = elem1[overlap_idx1, "Y"]         # Y from overlapping segment (same in both)
  
  eu = elem1[nonoverlap_idx1, "e"]      # e from non-overlapping segment of first element
  ev = elem2[nonoverlap_idx2, "e"]      # e from non-overlapping segment of second element
  ew = elem1[overlap_idx1, "e"]         # e from overlapping segment
  
  au = elem1[nonoverlap_idx1, "a"]      # a from non-overlapping segment of first element
  av = elem2[nonoverlap_idx2, "a"]      # a from non-overlapping segment of second element
  aw = elem1[overlap_idx1, "a"]         # a from overlapping segment
  
  r1 = elem1[nonoverlap_idx1, "r"]      # r from non-overlapping segment of first element
  r2 = elem1[overlap_idx1, "r"]         # r from overlapping segment in first element
  r3 = elem2[nonoverlap_idx2, "r"]      # r from non-overlapping segment of second element
  r4 = elem2[overlap_idx2, "r"]         # r from overlapping segment in second element
  
  # Derivative orders and evaluation point
  tu0 <- tv0 <- tw0 <- -oData$Segments$Time
  
  #--------------------------------
  # mgf-marginalisation computation
  #--------------------------------
  # Denominators at eval point
  LA <- beta  - r1*eu*tu0 - r2*ew*tw0
  LB <- beta  - r3*ev*tv0 - r4*ew*tw0
  Mu <- beta2 - au*tu0
  Mv <- beta2 - av*tv0
  Mw <- beta2 - aw*tw0
  
  # ================================
  # 1b. Vectorised log-scale coefficient method
  # ================================
  
  # Precompute log-factorials
  lfact <- lgamma(0:max(yu, yv, yw)+1)
  
  # m, n, p, q, s, t, u indices
  m <- 0:yu
  n <- 0:yw
  p <- 0:yv
  q <- 0:yw
  s <- 0:yu
  t <- 0:yv
  u <- 0:yw
  
  # Compute log(Atilde)
  M <- outer(m, n, "+")
  lfact_matA <- matrix(lfact[m+1], nrow=length(m), ncol=length(n)) +
    matrix(lfact[n+1], nrow=length(m), ncol=length(n), byrow=TRUE)
  logAtilde <- alpha*log(beta) + lgamma(alpha + M) - lgamma(alpha) +
    outer(m, n, function(mm, nn) mm*log(r1*eu) + nn*log(r2*ew)) -
    lfact_matA +
    (-alpha - M) * log(LA)
  
  # Compute log(Btilde)
  P <- outer(p, q, "+")
  lfact_matB <- matrix(lfact[p+1], nrow=length(p), ncol=length(q)) +
    matrix(lfact[q+1], nrow=length(p), ncol=length(q), byrow=TRUE)
  logBtilde <- alpha*log(beta) + lgamma(alpha + P) - lgamma(alpha) +
    outer(p, q, function(pp, qq) pp*log(r3*ev) + qq*log(r4*ew)) -
    lfact_matB +
    (-alpha - P) * log(LB)
  
  # Compute log(Cu, Cv, Cw)
  logCu <- alpha2*log(beta2) + lgamma(alpha2 + s) - lgamma(alpha2) + s*log(au) - (alpha2 + s)*log(Mu) - lfact[s+1]
  logCv <- alpha2*log(beta2) + lgamma(alpha2 + t) - lgamma(alpha2) + t*log(av) - (alpha2 + t)*log(Mv) - lfact[t+1]
  logCw <- alpha2*log(beta2) + lgamma(alpha2 + u) - lgamma(alpha2) + u*log(aw) - (alpha2 + u)*log(Mw) - lfact[u+1]
  
  # Vectorised accumulation in log-space
  logcoef_sum <- -Inf  # log(0) initialisation
  for (n_idx in 1:length(n)) {
    n_val <- n[n_idx]
    for (q_val in 0:(yw - n_val)) {
      sum_u <- logplusvec(logAtilde[, n_idx] + rev(logCu))  # sum over m
      sum_v <- logplusvec(logBtilde[, q_val+1] + rev(logCv)) # sum over p
      uw <- yw - n_val - q_val
      logcoef_sum <- logplus(logcoef_sum, sum_u + sum_v + logCw[uw+1])
    }
  }
  
  # Multiply by factorials in log-space
  log_deriv_coef <- logcoef_sum + lfact[yu+1] + lfact[yv+1] + lfact[yw+1]
  
  #---------------
  # CI computation
  #---------------
  A <- writeAmat(currentG)
  bkgdA <- bkgdAmat(currentG, iterations)
  Y <- extractYvec(currentG)
  log_p <- log_deriv_coef - sum(lfactorial(Y)) + sum(Y * log(oData$Segments$Time)) #!
  cat("p parameter: ", exp(log_p), "\n")
  sourceT_rates <- matrix(rgamma(Nsrc * iterations, shape = alpha, rate = beta / oData$Segments$Time), nrow = Nsrc)
  bkgdT_rates <- matrix(rgamma(Nseg * iterations, shape = bkgdShape, rate = bkgdRates / oData$Segments$Time), ncol = Nseg)
  bkgdATrates <- c(oData$Background$area * bkgdT_rates)
  region_obs <- matrix(rpois(Nseg * iterations, lambda = t(A %*% sourceT_rates) + bkgdA * bkgdT_rates), ncol = Nseg)
  bkgd_obs <- rpois(iterations, bkgdATrates)
  hist(region_obs[, 1], xlab = "simulated counts in segment 1", main = expression("Histogram of cluster ", binary_indices[i]))
  abline(v = Y[1], col = "red")
  hist(region_obs[, 2], xlab = "simulated counts in segment 2", main = expression("Histogram of cluster ", binary_indices[i]))
  abline(v = Y[2], col = "red")
  hist(region_obs[, 3], xlab = "simulated counts in segment 3", main = expression("Histogram of cluster ", binary_indices[i]))
  abline(v = Y[3], col = "red")
  region_y <- subset(region_obs, region_obs[, 1] == Y[1] & region_obs[, 2] == Y[2] & region_obs[, 3] == Y[3])
  cat("matching appearances: ", nrow(region_y), "\n") # number of times the observed Y happen
  cat("matching frequency", nrow(region_y) / iterations, "\n") # frequency the observed Y happens
  (freqq <- qbinom(c(0.025, 0.975), size = iterations, prob = exp(log_p))) # theoratical exact confidence interval
  cat("CI: ", freqq, "\n")
  inCI[i] <- nrow(region_y) >= freqq[1] & nrow(region_y) <= freqq[2] # the frequency Y happens is within the CI
  boundarySim[i] <- nrow(region_y) == 0
  boundaryCI[i] <- freqq[1] == 0
}

sum(inCI) / length(inCI)
sum(boundarySim) / length(boundarySim)
sum(boundaryCI) / length(boundaryCI)