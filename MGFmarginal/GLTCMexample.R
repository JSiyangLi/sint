# an example demonstrating the usefulness of the coefficient-Equating-convolution method.
# ================================
# Verification: coefficient method, Leibniz, caracas
# ================================

library(caracas)

poch <- function(a, n) {
  if (n == 0) return(1)
  gamma(a + n) / gamma(a)
}

# Parameters
beta  <- 3
beta2 <- 2
alpha <- 2
alpha2 <- 1

r1 <- 1; r2 <- 2; r3 <- 1; r4 <- 1
eu <- 1; ev <- 1; ew <- 1
au <- 1; av <- 1; aw <- 1

# Derivative orders and evaluation point
yu <- 2; yv <- 1; yw <- 2
tu0 <- tv0 <- tw0 <- -1

# Denominators at eval point
LA <- beta  - r1*eu*tu0 - r2*ew*tw0
LB <- beta  - r3*ev*tv0 - r4*ew*tw0
Mu <- beta2 - au*tu0
Mv <- beta2 - av*tv0
Mw <- beta2 - aw*tw0

# ================================
# 1. Coefficient method
# ================================

Atilde <- matrix(0, yu+1, yw+1)
for (m in 0:yu) {
  for (n in 0:yw) {
    Atilde[m+1, n+1] <- beta^alpha *
      poch(alpha, m + n) / (factorial(m) * factorial(n)) *
      (r1*eu)^m * (r2*ew)^n * LA^(-alpha - m - n)
  }
}

Btilde <- matrix(0, yv+1, yw+1)
for (p in 0:yv) {
  for (q in 0:yw) {
    Btilde[p+1, q+1] <- beta^alpha *
      poch(alpha, p + q) / (factorial(p) * factorial(q)) *
      (r3*ev)^p * (r4*ew)^q * LB^(-alpha - p - q)
  }
}

Cu <- sapply(0:yu, function(s)
  beta2^alpha2 * poch(alpha2, s) / factorial(s) * au^s * Mu^(-alpha2 - s)
)
Cv <- sapply(0:yv, function(t)
  beta2^alpha2 * poch(alpha2, t) / factorial(t) * av^t * Mv^(-alpha2 - t)
)
Cw <- sapply(0:yw, function(u)
  beta2^alpha2 * poch(alpha2, u) / factorial(u) * aw^u * Mw^(-alpha2 - u)
)

coef_sum <- 0
for (n in 0:yw) {
  for (q in 0:(yw - n)) {
    sum_u <- sum(sapply(0:yu, function(m) Atilde[m+1, n+1] * Cu[yu-m+1]))
    sum_v <- sum(sapply(0:yv, function(p) Btilde[p+1, q+1] * Cv[yv-p+1]))
    uw <- yw - n - q
    coef_sum <- coef_sum + sum_u * sum_v * Cw[uw+1]
  }
}
deriv_coef <- factorial(yu) * factorial(yv) * factorial(yw) * coef_sum

# ================================
# 2. Direct Leibniz (brute force, only small orders)
# ================================
direct_sum <- 0
for (m in 0:yu) {
  s <- yu - m
  for (p in 0:yv) {
    t <- yv - p
    for (n in 0:yw) {
      for (q in 0:(yw - n)) {
        u <- yw - n - q
        
        mult_u <- factorial(yu) / (factorial(m) * factorial(s))
        mult_v <- factorial(yv) / (factorial(p) * factorial(t))
        mult_w <- factorial(yw) / (factorial(n) * factorial(q) * factorial(u))
        multinom <- mult_u * mult_v * mult_w
        
        DA  <- beta^alpha * poch(alpha, m + n) * (r1*eu)^m * (r2*ew)^n * LA^(-alpha - m - n)
        DB  <- beta^alpha * poch(alpha, p + q) * (r3*ev)^p * (r4*ew)^q * LB^(-alpha - p - q)
        DCu <- beta2^alpha2 * poch(alpha2, s) * au^s * Mu^(-alpha2 - s)
        DCv <- beta2^alpha2 * poch(alpha2, t) * av^t * Mv^(-alpha2 - t)
        DCw <- beta2^alpha2 * poch(alpha2, u) * aw^u * Mw^(-alpha2 - u)
        
        direct_sum <- direct_sum + multinom * DA * DB * DCu * DCv * DCw
      }
    }
  }
}
deriv_leibniz <- direct_sum

# ================================
# 3. Symbolic with caracas (SymPy backend)
# ================================

tu <- symbol("tu"); tv <- symbol("tv"); tw <- symbol("tw")

f_sym <- (beta/(beta - r1*eu*tu - r2*ew*tw))^alpha *
  (beta/(beta - r3*ev*tv - r4*ew*tw))^alpha *
  (beta2/(beta2 - au*tu))^alpha2 *
  (beta2/(beta2 - av*tv))^alpha2 *
  (beta2/(beta2 - aw*tw))^alpha2

# Differentiate step by step (since der() is only first order)
df_sym <- f_sym
for (i in 1:yu) df_sym <- der(df_sym, tu)
for (i in 1:yv) df_sym <- der(df_sym, tv)
for (i in 1:yw) df_sym <- der(df_sym, tw)

#df_sym <- simplify(df_sym)

# Substitute numeric values and coerce to R numeric
df_eval <- subs(df_sym, list(tu = tu0, tv = tv0, tw = tw0))
deriv_caracas <- as_func(df_eval)()

# ================================
# 1b. Vectorised log-scale coefficient method
# ================================
library(Rcpp)
sourceCpp("logsum.cpp")

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
deriv_coef_log <- exp(log_deriv_coef)


# ================================
# 1e. Fully log-scale FFT-based coefficient method
# ================================

# Helper: next power of 2 for FFT padding
nextpow2 <- function(x) 2^ceiling(log2(x))

# 1. Compute log-coefficients for each factor
# Factor 1: Atilde (powers of tw: m+n)
m <- 0:yu; n <- 0:yw
M <- outer(m, n, "+")
lfact_matA <- matrix(lfact[m+1], nrow=length(m), ncol=length(n)) +
  matrix(lfact[n+1], nrow=length(m), ncol=length(n), byrow=TRUE)
logAtilde <- alpha*log(beta) + lgamma(alpha + M) - lgamma(alpha) +
  outer(m, n, function(mm, nn) mm*log(r1*eu) + nn*log(r2*ew)) -
  lfact_matA + (-alpha - M)*log(LA)
logA_tw <- sapply(1:(yw+1), function(j) logplusvec(logAtilde[,j] + rev(logCu)))

# Factor 2: Btilde (powers of tw: p+q)
p <- 0:yv; q <- 0:yw
P <- outer(p, q, "+")
lfact_matB <- matrix(lfact[p+1], nrow=length(p), ncol=length(q)) +
  matrix(lfact[q+1], nrow=length(p), ncol=length(q), byrow=TRUE)
logBtilde <- alpha*log(beta) + lgamma(alpha + P) - lgamma(alpha) +
  outer(p, q, function(pp, qq) pp*log(r3*ev) + qq*log(r4*ew)) -
  lfact_matB + (-alpha - P)*log(LB)
logB_tw <- sapply(1:(yw+1), function(j) logplusvec(logBtilde[,j] + rev(logCv)))

# Factor 3: Cw (powers of tw)
logCw_vec <- logCw  # already vector

# 2. FFT convolution in log-space
# Convert logs to normalized linear scale to avoid overflow
scale <- max(c(logA_tw, logB_tw, logCw_vec))
A_lin <- exp(logA_tw - scale)
B_lin <- exp(logB_tw - scale)
C_lin <- exp(logCw_vec - scale)

# Pad sequences for FFT
nfft <- nextpow2(length(A_lin) + length(B_lin) + length(C_lin) - 2)
A_fft <- fft(c(A_lin, rep(0, nfft - length(A_lin))))
B_fft <- fft(c(B_lin, rep(0, nfft - length(B_lin))))
C_fft <- fft(c(C_lin, rep(0, nfft - length(C_lin))))

# Multiply in frequency domain and inverse FFT
conv_fft <- Re(fft(A_fft * B_fft * C_fft, inverse = TRUE)) / nfft

# Pick coefficient of tw^yw (index yw+1)
coef_tw_yw <- conv_fft[yw+1]

# Rescale back to log-scale and multiply by factorials
deriv_fft <- coef_tw_yw * exp(3*scale) * factorial(yu) * factorial(yv) * factorial(yw)


# ================================
# Compare
# ================================
cat("Derivative (coefficient method):", deriv_coef, "\n")
cat("Derivative (direct Leibniz)    :", deriv_leibniz, "\n")
cat("Derivative (caracas symbolic)  :", deriv_caracas, "\n")
cat("Derivative (coefficient method, log-scale):", deriv_coef_log, "\n")
cat("Derivative (coefficient method, log-scale FFT):", deriv_fft, "\n")
cat("all.equal (coef, Leibniz, caracas, log-scale coef, log-scale coef FFT)    :", 
    all.equal(deriv_coef, deriv_leibniz, deriv_caracas, deriv_coef_log, deriv_fft), "\n")
