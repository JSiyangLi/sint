library(Rcpp)
library(kStatistics)
library(tictoc)
sourceCpp("logsum.cpp")
#sourceCpp("combo_iterator.cpp")

#########################
# hierarchical bkgd model
#########################
# for isolated sources
mgfD <- function( # for each single source
    y, # the number of photons observed
    alpha, # source intensity population alpha
    beta, # source intensity population beta
    e, # exposure of the isolated source
    a, # area of the isolates source
    Time, # observation time
    bkgdalpha, # local background population alpha
    bkgdbeta # local background population beta
) {
  ##################
  # helper functions
  ##################
  
  #++++++++++++++++
  # cgf derivatives
  #++++++++++++++++
  cgf_derivative_individual_log <- function(
    alpha, beta, k, e, a, Time, bkgdalpha, bkgdbeta 
  ) {
    # Compute first term: -alpha * (k-1)! * (beta - 0.9*e*Time)^(-k) * (0.9*e)^k
    term1_inner <- logplus(log(beta), log(0.9) + log(e) + log(Time))
    term1_log <- log(alpha) + lfactorial(k - 1) - k * term1_inner + k * (log(0.9) + log(e))
    
    # Compute second term: -bkgdalpha * (k-1)! * (bkgdbeta - a*Time)^(-k) * a^k
    term2_inner <- logplus(log(bkgdbeta), log(a) + log(Time))
    term2_log <- log(bkgdalpha) + lfactorial(k - 1) - k * term2_inner + k * log(a)
    
    # Combine terms (both always negative)
    log_abs_result <- logplus(term1_log, term2_log)
    
    return(log_abs_result)  # Always positive
  }
  
  #++++++++++++++++
  # Bell polynomial
  #++++++++++++++++
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
  
  #++++
  # mgf
  #++++
  lmgf <- function(alpha, beta, e, a, Time, bkgdalpha, bkgdbeta) {
    term1 <- alpha * log(beta) - alpha * logplus(log(beta), log(0.9) + log(e) + log(Time))
    term2 <- bkgdalpha * log(bkgdbeta) - bkgdalpha * logplus(log(bkgdbeta), log(a) + log(Time))
    term1 + term2
  }
  
  #############
  # computation
  #############
  if (y == 0) {
    l0mgf <- lmgf(alpha=alpha, beta=beta, e=e, a=a, Time=Time, bkgdalpha=bkgdalpha, bkgdbeta=bkgdbeta)
    return(list(logDeriv = l0mgf, sign = 1))
  } else if (y >= 1) {
    lDcgf <- sapply(1:y, 
                    FUN = function(i) {cgf_derivative_individual_log(alpha=alpha, beta=beta, k=i, e=e, a=a, Time=Time, bkgdalpha=bkgdalpha, bkgdbeta=bkgdbeta)})
    lDBell <- recursive_logBellPol_signed(n=y, logv = lDcgf, vsign = rep(1, y)) # the cgf derivatives are always positive
    l0mgf <- lmgf(alpha=alpha, beta=beta, e=e, a=a, Time=Time, bkgdalpha=bkgdalpha, bkgdbeta=bkgdbeta)
    return(list(logDeriv = l0mgf + lDBell$logB, sign = lDBell$sign))
  } else { # y < 0
    stop("invalid Poisson observation")
  }
  
}

lmgfMarg <- function(
    SData, # Separated iSolated Source info in the given format
    alpha, # source intensity population alpha
    beta, # source intensity population beta
    Time, # observation time
    bkgdalpha, # local background population alpha
    bkgdbeta # local background population beta
    ) {
  yea <- cbind(SData$Isolate$cnt.obs, SData$Isolate$eff.rea, SData$Isolate$area) # matrix of y, e and a
  Sresult <- apply(yea, 1, FUN = function(v) {mgfD(y=v[1], alpha, beta, e=v[2], a=v[3], Time, bkgdalpha, bkgdbeta)}) |> unlist() # iterating mgfD
  even_indices <- 2 * (1:(length(Sresult)/2)) # the result is written in a vector where even indices are signs and odd indices are logged absolute values
  
  # extract the logged values and signs, and compute the density
  lProdsign <- prod(Sresult[even_indices])
  if (lProdsign <= 0) { # densities should always be positive
    warning("Implementation error: mgf-derived density computed to be negative.")
  }
  lSumderiv <- sum(Sresult[even_indices-1]) - sum(lfactorial(yea[, 1])) + sum(yea[, 1] * log(Time)) # mgf marginalisation equation
  
  return(lSumderiv)
}




#####################
# overlapping sources
#####################
log_rising_factorial <- function(x, n) {
  if (n == 0) return(0)
  lgamma(x + n) - lgamma(x)
}


mixed_derivative_log <- function(alpha, beta, yu, yv, yw,
                                 r1, eu, r2, ew,
                                 r3, ev, r4,
                                 alpha2, beta2, au, av, aw,
                                 tu0, tv0, tw0) {
  # Helper functions
  log_rising_factorial <- function(x, n) {
    if (n == 0) return(0)
    lgamma(x + n) - lgamma(x)
  }
  
  # Precompute denominators
  denom_f1 <- beta - tu0 * r1 * eu - tw0 * r2 * ew
  denom_f2 <- beta - tv0 * r3 * ev - tw0 * r4 * ew
  denom_g3 <- beta2 - aw * tw0
  
  # Initialize total log-abs and sign
  total_log_abs <- -Inf
  total_sign <- 1
  
  # Main loop over i (tw derivatives for f1) and j (tw derivatives for f2)
  for (i in 0:yw) {
    for (j in 0:(yw - i)) {
      k3 <- yw - i - j
      
      # Multinomial coefficient
      log_multinomial <- lfactorial(yw) - (lfactorial(i) + lfactorial(j) + lfactorial(k3))
      
      # g3 component
      log_g3 <- k3 * log(aw) + 
        log_rising_factorial(alpha2, k3) + 
        (alpha2 + k3) * (log(beta2) - log(denom_g3))
      
      # Term1: f1 and g1 components
      term1_log_abs <- -Inf
      term1_sign <- 1
      for (k in 0:yu) {
        n1 <- yu - k
        log_term <- (
          lchoose(yu, k) +
            k * log(au) +
            log_rising_factorial(alpha2, k) +
            n1 * log(r1 * eu) +
            i * log(r2 * ew) +
            log_rising_factorial(alpha, n1 + i) +
            alpha * log(beta) - 
            (alpha + n1 + i) * log(denom_f1)
        )
        # Update term1 using logplus
        if (is.infinite(term1_log_abs)) {
          term1_log_abs <- log_term
        } else {
          sum_log <- logplus(term1_log_abs, log_term)
          term1_log_abs <- sum_log
        }
      }
      
      # Term2: f2 and g2 components
      term2_log_abs <- -Inf
      term2_sign <- 1
      for (l in 0:yv) {
        n2 <- yv - l
        log_term <- (
          lchoose(yv, l) +
            l * log(av) +
            log_rising_factorial(alpha2, l) +
            n2 * log(r3 * ev) +
            j * log(r4 * ew) +
            log_rising_factorial(alpha, n2 + j) +
            alpha * log(beta) - 
            (alpha + n2 + j) * log(denom_f2)
        )
        # Update term2 using logplus
        if (is.infinite(term2_log_abs)) {
          term2_log_abs <- log_term
        } else {
          sum_log <- logplus(term2_log_abs, log_term)
          term2_log_abs <- sum_log
        }
      }
      
      # Combine all components
      current_log_abs <- log_multinomial + log_g3 + term1_log_abs + term2_log_abs
      
      # Update total using logplus (all terms positive)
      if (is.infinite(total_log_abs)) {
        total_log_abs <- current_log_abs
      } else {
        total_log_abs <- logplus(total_log_abs, current_log_abs)
      }
    }
  }
  
  return(list(log_abs = total_log_abs, sign = 1))
}
