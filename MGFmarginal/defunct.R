### from mgf.R
#######################
### mgf derivatives ###
#######################
# matrix-apply-convolve implementation
log_gamma_mgf_uni <- function(
    alpha, beta, t
) {
  if (any(beta <= 0)) {
    cat("invalid shape or rate values: alpha: ", alpha, "beta: ", beta)
    return(-Inf)
  }
  
  alpha * c(log(beta) - log(beta - t)) # recycling the betas for different alphas
}

log_gamma_bkgd_part <- function( # background part in the overall derivative
  umax,
  v,
  k, # the index of which segments are involved (single-source)
  response, # response vector
  design, # design matrix
  Time # observation time
) {
  alpha.xi <- 1e-06 # values from Lazhi's paper
  n <- nrow(design)
  Il <- ncol(design) - 1
  mgf.t <- rep(-Time, n) # mgf input vector
  beta.xi <- Time * design[n, Il + 1] * 1e-12
  x <- response[length(response)]
  univariate.t <- t(mgf.t) %*% design[, ncol(design)]
  
  usum <- (0:umax) + v # summands in the convolution
  a <- design[k, ncol(design)] # vector of areas for each segment
  
  - usum * beta.xi + lgamma(alpha.xi + x + usum) - lgamma(alpha.xi + x) + log_gamma_mgf_uni(alpha.xi + x + usum, beta.xi, univariate.t)
}


log_gamma_bkgd <- function(
    v=0, # added constant from the secondary source
    k, # the index of which segments are involved (single-source). if length == 1, no secondary segment considered
    response, # response vector
    design, # design matrix
    Time # observation time
) {
  if (length(v) == 1) {
    if (v != 0) {
      if (length(k) != length(v) + 1) {
        cat("please check segment indices")
        return(NaN)
      }
    }
  }
  
  umax <- response[k[1]]
  usum <- 0:umax # summands in the convolution
  
  log_gamma_bkgd_front <- function(k, u_single) {
    a <- design[k, ncol(design)] # vector of areas for each segment
    u_single * log(a)
  }
  front <- ifelse(
    length(k) == 1, 
    log_gamma_bkgd_front(k[1], usum), 
    log_gamma_bkgd_front(k[1], usum) + sum(sapply(1:length(v), FUN = function(i) log_gamma_bkgd_front(k[i+1], v[i]))))
  front + log_gamma_bkgd_part(umax, sum(v), k[1], response, design, Time)
}

log_gamma_summand_1source <- function(
    k, # the index integer of which segments are involved (single-source)
    alpha, # population shape
    beta, # population rate
    design, # design matrix
    response,  # response vector
    Time # observation time
) {
  y <- response[k]
  u <- y:0 # !!! REVERSED order for convolution, representing y-u
  non0ind <- which(design[k, ] != 0) # non-0 indicies, including w and a, for single-source segment
  n <- nrow(design)
  r <- design[k, non0ind[1]]
  mgf.t <- rep(-Time, n) # mgf input vector
  univariate.t <- t(mgf.t) %*% design[, non0ind[1]]
  
  if (any(u < 0) | any(u > y)) {
    cat("invalid summand: ", u, " for observation: ", y)
    return(-Inf)
  }
  
  lchoose(y, u) + u * (log(r) - log(beta)) + lgamma(alpha + u) - lgamma(alpha) + log_gamma_mgf_uni(alpha + u, beta, univariate.t)
}
log_gamma_summand_1source(k=3, alpha = (1.2e-6)^2 / 8.7e-12, beta = 1.2e-6/8.7e-12, 
                          design = oDesignRes$Design, response = oDesignRes$Response, Time = oDesignRes$Time)

sapply(1:3, FUN = function(k) log_gamma_summand_1source(k, alpha = (1.2e-6)^2 / 8.7e-12, beta = 1.2e-6/8.7e-12, 
                                                        design = oDesignRes$Design, response = oDesignRes$Response, Time = oDesignRes$Time))

log_gamma_bkgd(v=1:2, k=1:3, response = xDesignRes$Response, design = xDesignRes$Design, Time = xDesignRes$Time)

sapply(0:34, FUN = function(w) {log_gamma_bkgd(v=w, k=1:2, response = xDesignRes$Response, design = xDesignRes$Design, Time = xDesignRes$Time)})

############
# convolutions of mgfs
###########
# expanding the count vector into a matrix
generate_count_matrix <- function(vvec) {
  if (length(vvec) > 3) {
    warning("the resulting matrix might be too large")
  }
  
  # Generate all possible combinations of counts
  grid <- expand.grid(lapply(vvec, function(x) 0:x))
  
  # Convert to matrix (and ensure integers)
  mat <- as.matrix(grid)
  storage.mode(mat) <- "integer"  # Ensure integer storage
  
  # Set column names (optional)
  colnames(mat) <- NULL  # Remove column names for cleaner output
  
  return(mat)
}

log_gamma_convolve1 <- function(
    k, # the index integer of which segments are involved (single-source)
    alpha = (1.2e-6)^2 / 8.7e-12, # population shape
    beta = 1.2e-6/8.7e-12, # population rate
    design = oDesignRes$Design, # design matrix
    response = oDesignRes$Response,  # response vector
    Time = oDesignRes$Time # observation time
) {
  if (length(k) != 1) {
    cat("Incorrect k length, use log_gamma_convolve2")
    return(NaN)
  }
  MGFsource <- log_gamma_summand_1source(k, alpha, beta, design, response, Time)
  MGFbkgd <- log_gamma_bkgd(v=0, k, response, design, Time)
  if (length(MGFsource) != length(MGFbkgd)) {
    cat("Incorrect MGF lengths")
    return(NaN)
  }
  
  MGFsource + MGFbkgd
}

log_gamma_convolve2 <- function(
    k, # the index integer of which segments are involved (single-source)
    alpha = (1.2e-6)^2 / 8.7e-12, # population shape
    beta = 1.2e-6/8.7e-12, # population rate
    design = xDesignRes$Design, # design matrix
    response = xDesignRes$Response,  # response vector
    Time = xDesignRes$Time # observation time
) {
  if (length(k) != 2) {
    cat("Incorrect k length, use log_gamma_convolve1")
    return(NaN)
  }
  MGFsource1 <- log_gamma_summand_1source(k[1], alpha, beta, design, response, Time)
  MGFsource2 <- log_gamma_summand_1source(k[2], alpha, beta, design, response, Time)
  vmax <- response[k[2]]
  vsum <- 0:vmax
  MGFbkgd <- sapply(vsum, FUN = function(w) {log_gamma_bkgd(v=w, k, response, design, Time)})
  if (length(MGFsource1) != nrow(MGFbkgd) | length(MGFsource2) != ncol(MGFbkgd)) {
    cat("Incorrect MGF lengths")
    return(NaN)
  }
  
  MGFsource1mat <- rep(MGFsource1, length(MGFsource2)) |> matrix(ncol = length(MGFsource2))
  ((MGFsource1mat + MGFbkgd) |> apply(MARGIN = 2, FUN = logplusvec)) + MGFsource2
}

log_gamma_convolve_array <- function(
    k, # the index integer of which segments are involved (single-source)
    alpha = (1.2e-6)^2 / 8.7e-12, # population shape
    beta = 1.2e-6/8.7e-12, # population rate
    design = xDesignRes$Design, # design matrix
    response = xDesignRes$Response,  # response vector
    Time = xDesignRes$Time # observation time
) {
  if (length(k) == 1) {
    cat("Incorrect k length, use log_gamma_convolve1")
    return(NaN)
  }
  
  MGFsource <- sapply(k, FUN = log_gamma_summand_1source, 
                      alpha=alpha, beta=beta, design=design, response=response, Time=Time)
  
  for (ki in 1:length(k)) {
    vmax <- response[ki]
    vsum <- 0:v
    MGFbkgd <- sapply(vsum, FUN = function(w) {log_gamma_bkgd(v=w, k, response, design, Time)})
  }
  for (v in vmax) {
    
    if (length(MGFsource1) != nrow(MGFbkgd) | length(MGFsource2) != ncol(MGFbkgd)) {
      cat("Incorrect MGF lengths")
      return(NaN)
    }
  }
  MGFsource1mat <- rep(MGFsource1, length(MGFsource2)) |> matrix(ncol = length(MGFsource2))
  ((MGFsource1mat + MGFbkgd) |> apply(MARGIN = 2, FUN = logplusvec)) + MGFsource2
}

log_gamma_convolve2(k=1:2)

log_gamma_recursion <- function(
    kmax, # the index integer of which segments are involved (single-source)
    alpha = (1.2e-6)^2 / 8.7e-12, # population shape
    beta = 1.2e-6/8.7e-12, # population rate
    design = oDesignRes$Design, # design matrix
    response = oDesignRes$Response,  # response vector
    Time = oDesignRes$Time # observation time
) {
  k = 2
  MGFsum = 0
  while (k < kmax) {
    y <- response[k]
    logsum_coef <- sapply(0:y, FUN = function(w) {log_gamma_convolve1(v=w, k-1)}) |> apply(MARGIN = 2, FUN = logplusvec)
    MGFsource <- log_gamma_summand_1source(k, alpha, beta, design, response, Time)
    MGFsum <- logsum_coef + MGFsource + MGFsum
    
    k = k + 1
  }
  
}


logplusvec(log_gamma_convolve1(v=10, k=1))
logsum_coef <- (sapply(0:34, FUN = function(w) {log_gamma_convolve(v=w, k=1:2)}) |> apply(MARGIN = 2, FUN = logplusvec)) + log_gamma_summand_1source(k=2, alpha = (1.2e-6)^2 / 8.7e-12, beta = 1.2e-6/8.7e-12, 
                                                                                                                                                     design = oDesignRes$Design, response = oDesignRes$Response, Time = oDesignRes$Time)


logsum_coef(1) + log_gamma_summand_1source(k=3, alpha = (1.2e-6)^2 / 8.7e-12, beta = 1.2e-6/8.7e-12, 
                                           design = oDesignRes$Design, response = oDesignRes$Response, Time = oDesignRes$Time)

convolve(log_gamma_summand_1source(k=1, alpha = (1.2e-6)^2 / 8.7e-12, beta = 1.2e-6/8.7e-12, 
                                   design = oDesignRes$Design, response = oDesignRes$Response, Time = oDesignRes$Time) |> exp(),
         log_gamma_bkgd(v=100, k=1:2, response = oDesignRes$Response, design = oDesignRes$Design, Time = oDesignRes$Time) |> exp(),
         type = "filter")




########################
# recursive implementation
########################
# Log-Gamma difference: log(Gamma(x + u) - log(Gamma(x)) = lgamma(x + u) - lgamma(x)
log_gamma_diff <- function(x, u) {
  lgamma(x + u) - lgamma(x)
}

log_gamma_summand_coef <- function(
    u,
    k, # the index integer of which segments are involved (single-source)
    alpha, # population shape
    beta, # population rate
    design, # design matrix
    response,  # response vector
    Time # observation time
) {
  y <- response[k]
  non0ind <- which(design[k, ] != 0) # non-0 indicies, including w and a, for single-source segment
  n <- nrow(design)
  r <- design[k, non0ind[1]]
  Il <- ncol(design) - 1
  mgf.t <- rep(-Time, n) # mgf input vector
  univariate.t <- t(mgf.t) %*% design[, non0ind[1]]
  bkgd.t <- t(mgf.t) %*% design[, ncol(design)]
  a <- design[k, ncol(design)] # vector of areas for each segment
  alpha.xi <- 1e-06 # values from Lazhi's paper
  beta.xi <- Time * design[n, Il + 1] * 1e-12
  x <- response[length(response)]
  
  if (any(u < 0) | any(u > y)) {
    cat("invalid summand: ", u, " for observation: ", y)
    return(-Inf)
  }
  
  lchoose(y, y-u) + (y-u) * (log(r) - log(beta)) + lgamma(alpha + y-u) - lgamma(alpha) + log_gamma_mgf_uni(alpha + y-u, beta, univariate.t) + u * (log(a) - log(beta.xi)) + log_gamma_mgf_uni(u, beta.xi, bkgd.t)
}
#log_gamma_summand_coef(u=1, k=3, alpha = (1.2e-6)^2 / 8.7e-12, beta = 1.2e-6/8.7e-12, 
#                          design = oDesignRes$Design, response = oDesignRes$Response, Time = oDesignRes$Time)


###################
# recursion computation
###################
compute_recursive_sum <- function(x, max_vals, g, memo = list()) {
  n <- length(max_vals)
  if (n == 0) return(0)  # Base case
  
  # Check memoization
  key <- paste(x, n, sep = ",")
  if (!is.null(memo[[key]])) return(memo[[key]])
  
  # Initialize log-sum for current layer
  log_total <- -Inf  # log(0)
  
  # Iterate over possible u_n values (0 to max_vals[n])
  for (u in 0:max_vals[n]) {
    # Compute term: log(g(u)) + log_gamma_diff(x, u)
    log_term <- g(u) + log_gamma_diff(x, u)
    
    # Recursive call for next layer (x + u)
    if (n > 1) {
      log_term <- log_term + compute_recursive_sum(x + u, max_vals[-n], g, memo)
    }
    
    # Accumulate in log-space
    log_total <- logplusvec(c(log_total, log_term))
  }
  
  # Memoize and return
  memo[[key]] <- log_total
  return(log_total)
}

result_log <- compute_recursive_sum(x = oData$Background$cnt + 1e-06, oData$Segments$cnt.obs[1:3], g = function(u) log_gamma_summand_coef(u, k=1, alpha = (1.2e-6)^2 / 8.7e-12, beta = 1.2e-6/8.7e-12, 
                                                                                                                                          design = oDesignRes$Design, response = oDesignRes$Response, Time = oDesignRes$Time))
# NegBin verification
lmgfMarg(xSData, alpha = (1.2e-6)^2 / 8.7e-12, beta = 1.2e-6/8.7e-12, Time = xData$Segments$Time, 
         bkgdalpha = 1e-06 + xData$Background$cnt, bkgdbeta = xData$Segments$Time * xData$Background$area * (1+1e-12)) - sum(lfactorial(xSData$Isolate$cnt.obs)) + sum(xSData$Isolate$cnt.obs * log(xData$Segments$Time))
mgfD(xSData$Isolate$cnt.obs[1], alpha = (1.2e-6)^2 / 8.7e-12, beta = 1.2e-6/8.7e-12, 
     e = xSData$Isolate$eff.rea[1], a = xSData$Isolate$area[1], Time = xData$Segments$Time, bkgdalpha = 1e-06 + xData$Background$cnt, bkgdbeta = xData$Segments$Time * xData$Background$area * (1+1e-12)) - lfactorial(xSData$Isolate$cnt.obs[1])
dnbinom(xSData$Isolate$cnt.obs[1], size = (1.2e-6)^2 / 8.7e-12 + 1e-06 + xData$Background$cnt, prob = (1.2e-6/8.7e-12) / (1 + 1.2e-6/8.7e-12), log = TRUE)

xsources <- dnbinom(0:xSData$Isolate$cnt.obs[1], (1.2e-6)^2 / 8.7e-12, sourceProbs[1])
xbkgds <- dnbinom(xSData$Isolate$cnt.obs[1]:0, 1e-06 + xData$Background$cnt, bkgdProbs[1])
convolve(xsources, xbkgds, type = "filter")

# "pure source"
dnbinom(xSData$Isolate$cnt.obs, size = sourceShape, prob = sourceProbs, log = TRUE) |> sum()
log_gamma_mgf_derivative <- function(order, alpha, beta, s) {
  lgamma(alpha + order) - lgamma(alpha) + alpha * log(beta) - (alpha + order) * logplus(log(beta),  log(s))
}
log_gamma_mgf_derivative(xSData$Isolate$cnt.obs[1], alpha = sourceShape, beta = sourceRates, s = xData$Segments$Time * 0.9 * xSData$Isolate$eff.rea[1])
log_marginalisation <- function(alpha, beta, erT, y) {
  sum(y * log(erT) - lfactorial(y) + sapply(1:length(y), function(i) log_gamma_mgf_derivative(order = y[i], alpha, beta, erT[i])))
}
log_marginalisation(y = xSData$Isolate$cnt.obs, alpha = sourceShape, beta = sourceRates, erT = xData$Segments$Time * 0.9 * xSData$Isolate$eff.rea)

log_gamma_mgf_derivative2 <- function(order, alpha, beta, s, r) {
  lgamma(alpha + order) - lgamma(alpha) + order * log(r) + alpha * log(beta) - (alpha + order) * logplus(log(beta),  log(s) + log(r))
}
log_marginalisation2 <- function(alpha, beta, Time, er, y) {
  sum(y * log(Time) - lfactorial(y) + sapply(1:length(y), function(i) log_gamma_mgf_derivative2(order = y[i], alpha, beta, s = Time, r = er[i])))
}
log_marginalisation2(y = xSData$Isolate$cnt.obs, alpha = sourceShape, beta = sourceRates, Time = xData$Segments$Time, er = 0.9 * xSData$Isolate$eff.rea)

# "pure background"
dnbinom(xSData$Isolate$cnt.obs, size = bkgdShape, prob = bkgdProbs, log = TRUE) |> sum()
#b1 <- sapply(1:length(xSData$Isolate$cnt.obs), function(i) log_gamma_mgf_derivative(order = xSData$Isolate$cnt.obs[i], alpha = bkgdShape, beta = bkgdRates, xData$Segments$Time * xSData$Isolate$area[i]))
#b2 <- log(xData$Segments$Time * xSData$Isolate$area) * xSData$Isolate$cnt.obs
#sum(b2 - lfactorial(xSData$Isolate$cnt.obs) + b1)
log_marginalisation(y = xSData$Isolate$cnt.obs, alpha = bkgdShape, beta = bkgdRates, erT = xData$Segments$Time * xSData$Isolate$area)

# Bell polynomial backups
log_gamma_mgf_derivative_bell2 <- function(order, alpha, beta, s, r, bkgdalpha, bkgdbeta, a) {
  # Base case: 0th derivative is the MGF itself
  if (order == 0) {
    logM0 <- alpha * log(beta) - alpha * logplus(log(beta), log(s) + log(r)) + bkgdalpha * log(bkgdbeta) - bkgdalpha * logplus(log(bkgdbeta), log(s) + log(a))
    return(list(logDeriv = logM0, sign = 1))
  }
  
  recursive_logBellPol_signed <- function(n, logv = rep(0, n), vsign = rep(1, n)) {
    # Input validation
    stopifnot(length(logv) >= n, length(vsign) >= n)
    
    # Initialize logB (log-abs values) and signB (signs)
    logB <- c(0, rep(-Inf, n))  # logB[1] = log(1) = 0 (B₀ = 1)
    signB <- c(1, rep(1, n))     # signB[1] = 1
    
    for (i in 1:n) {
      pos_terms <- neg_terms <- rep(-Inf, i)
      pos_idx <- neg_idx <- 1
      
      for (k in 1:i) {
        # Compute log-term: log|C(i-1,k-1) * v_k * B_{i-k}|
        log_term <- lchoose(i - 1, k - 1) + logv[k] + logB[i - k + 1]
        term_sign <- vsign[k] * signB[i - k + 1]
        
        if (term_sign > 0) {
          pos_terms[pos_idx] <- log_term
          pos_idx <- pos_idx + 1
        } else if (term_sign < 0) {
          neg_terms[neg_idx] <- log_term
          neg_idx <- neg_idx + 1
        }
        # term_sign == 0: Skip (log_term = -Inf)
      }
      
      # Sum positive/negative terms (handle empty cases)
      sum_pos <- if (pos_idx > 1) logplusvec(pos_terms[1:(pos_idx - 1)]) else -Inf
      sum_neg <- if (neg_idx > 1) logplusvec(neg_terms[1:(neg_idx - 1)]) else -Inf
      
      # Determine logB and signB (fixed logic)
      if (is.infinite(sum_pos) && is.infinite(sum_neg)) {
        logB[i + 1] <- -Inf  # All terms zero
        signB[i + 1] <- 1
      } else if (is.infinite(sum_neg)) {
        logB[i + 1] <- sum_pos  # Only positive terms
        signB[i + 1] <- 1
      } else if (is.infinite(sum_pos)) {
        logB[i + 1] <- sum_neg  # Only negative terms
        signB[i + 1] <- -1
      } else {
        # Mixed terms: Use logminus
        if (sum_pos >= sum_neg) {
          logB[i + 1] <- logminus(sum_pos, sum_neg)
          signB[i + 1] <- 1
        } else {
          logB[i + 1] <- logminus(sum_neg, sum_pos)
          signB[i + 1] <- -1
        }
      }
    }
    
    return(list(logB = logB[n + 1], sign = signB[n + 1]))
  }
  
  # Helper function for Gamma CGF derivatives
  gamma_cgf_deriv_log <- function(k, alpha, beta, s, r, bkgdalpha, bkgdbeta, a) {
    # k-th derivative of CGF at s: K^{(k)}(s) = alpha * (k-1)! / (beta - s)^k
    logsource <- log(alpha) + k * log(r) + lfactorial(k - 1) - k * logplus(log(beta), log(s) + log(r))
    logbkgd <- log(bkgdalpha) + k * log(a) + lfactorial(k - 1) - k * logplus(log(bkgdbeta), log(s) + log(a))
    logplus(logsource, logbkgd)
  }
  
  # Compute CGF derivatives (orders 1 to 'order')
  logv <- sapply(1:order, function(k) gamma_cgf_deriv_log(k, alpha, beta, s, r, bkgdalpha, bkgdbeta, a))
  vsign <- rep(1, order)  # All derivatives are positive
  
  # Compute Bell polynomial for the given order
  bell_res <- recursive_logBellPol_signed(n = order, logv = logv, vsign = vsign)
  cat("bell_res", bell_res$logB, "\n")
  # Compute log(M(s)) = K(s)
  logM <- alpha * log(beta) - alpha * logplus(log(beta), log(s) + log(r)) + bkgdalpha * log(bkgdbeta) - bkgdalpha * logplus(log(bkgdbeta), log(s) + log(a))
  cat("mgf", logM, "\n")
  # Final derivative: M^{(n)}(s) = M(s) * B_n
  list(logDeriv = logM + bell_res$logB, sign = bell_res$sign)
}

log_gamma_mgf_derivative_bell2(xSData$Isolate$cnt.obs[1], alpha = sourceShape, beta = sourceRates, s = xData$Segments$Time, r = 0.9 * xSData$Isolate$eff.rea[1],
                               bkgdalpha = bkgdShape, bkgdbeta = bkgdRates, a = xSData$Isolate$area[1])
# Updated marginalization function using Bell polynomial approach
log_marginalisation_bell2 <- function(alpha, beta, er, Time, y, bkgdalpha, bkgdbeta, a) {
  term_sum <- sum(y * log(Time) - lfactorial(y))
  
  bell_terms <- sapply(seq_along(y), function(i) {
    deriv_res <- log_gamma_mgf_derivative_bell2(
      order = y[i], 
      alpha = alpha,
      beta = beta,
      s = Time,  # Evaluate at negative exposure
      r = er[i],
      bkgdalpha = bkgdalpha,
      bkgdbeta = bkgdbeta,
      a = a[i]
    )
    # Safety check (derivatives should be positive)
    if (deriv_res$sign < 0) {
      warning("Unexpected negative derivative sign at i=", i)
    }
    deriv_res$logDeriv
  })
  
  term_sum + sum(bell_terms)
}
log_marginalisation_bell2(y = xSData$Isolate$cnt.obs[1], alpha = sourceShape, beta = sourceRates, Time = xData$Segments$Time, er = 0.9 * xSData$Isolate$eff.rea[1],
                          bkgdalpha = bkgdShape, bkgdbeta = bkgdRates, a = xSData$Isolate$area[1])

#!!!!!!!!!!functional backup
log_gamma_mgf_derivative_bell <- function(order, alpha, beta, s) {
  # Base case: 0th derivative is the MGF itself
  if (order == 0) {
    logM0 <- alpha * log(beta) - alpha * logplus(log(beta), log(s))
    return(list(logDeriv = logM0, sign = 1))
  }
  
  recursive_logBellPol_signed <- function(n, logv = rep(0, n), vsign = rep(1, n)) {
    # Input validation
    stopifnot(length(logv) >= n, length(vsign) >= n)
    
    # Initialize logB (log-abs values) and signB (signs)
    logB <- c(0, rep(-Inf, n))  # logB[1] = log(1) = 0 (B₀ = 1)
    signB <- c(1, rep(1, n))     # signB[1] = 1
    
    for (i in 1:n) {
      pos_terms <- neg_terms <- rep(-Inf, i)
      pos_idx <- neg_idx <- 1
      
      for (k in 1:i) {
        # Compute log-term: log|C(i-1,k-1) * v_k * B_{i-k}|
        log_term <- lchoose(i - 1, k - 1) + logv[k] + logB[i - k + 1]
        term_sign <- vsign[k] * signB[i - k + 1]
        
        if (term_sign > 0) {
          pos_terms[pos_idx] <- log_term
          pos_idx <- pos_idx + 1
        } else if (term_sign < 0) {
          neg_terms[neg_idx] <- log_term
          neg_idx <- neg_idx + 1
        }
        # term_sign == 0: Skip (log_term = -Inf)
      }
      
      # Sum positive/negative terms (handle empty cases)
      sum_pos <- if (pos_idx > 1) logplusvec(pos_terms[1:(pos_idx - 1)]) else -Inf
      sum_neg <- if (neg_idx > 1) logplusvec(neg_terms[1:(neg_idx - 1)]) else -Inf
      
      # Determine logB and signB (fixed logic)
      if (is.infinite(sum_pos) && is.infinite(sum_neg)) {
        logB[i + 1] <- -Inf  # All terms zero
        signB[i + 1] <- 1
      } else if (is.infinite(sum_neg)) {
        logB[i + 1] <- sum_pos  # Only positive terms
        signB[i + 1] <- 1
      } else if (is.infinite(sum_pos)) {
        logB[i + 1] <- sum_neg  # Only negative terms
        signB[i + 1] <- -1
      } else {
        # Mixed terms: Use logminus
        if (sum_pos >= sum_neg) {
          logB[i + 1] <- logminus(sum_pos, sum_neg)
          signB[i + 1] <- 1
        } else {
          logB[i + 1] <- logminus(sum_neg, sum_pos)
          signB[i + 1] <- -1
        }
      }
    }
    
    return(list(logB = logB[n + 1], sign = signB[n + 1]))
  }
  
  # Helper function for Gamma CGF derivatives
  gamma_cgf_deriv_log <- function(k, alpha, beta, s) {
    # k-th derivative of CGF at s: K^{(k)}(s) = alpha * (k-1)! / (beta - s)^k
    log(alpha) + lfactorial(k - 1) - k * logplus(log(beta), log(s))
  }
  
  # Compute CGF derivatives (orders 1 to 'order')
  logv <- sapply(1:order, function(k) gamma_cgf_deriv_log(k, alpha, beta, s))
  vsign <- rep(1, order)  # All derivatives are positive
  
  # Compute Bell polynomial for the given order
  bell_res <- recursive_logBellPol_signed(n = order, logv = logv, vsign = vsign)
  
  # Compute log(M(s)) = K(s)
  logM <- alpha * log(beta) - alpha * logplus(log(beta), log(s))
  
  # Final derivative: M^{(n)}(s) = M(s) * B_n
  list(logDeriv = logM + bell_res$logB, sign = bell_res$sign)
}
