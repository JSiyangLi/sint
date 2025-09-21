library(fftwtools)

# Your fftn() implementation
fftn_working <- function(x, inverse = FALSE) {
  dims <- dim(x)
  X <- x
  for (axis in seq_along(dims)) {
    X <- apply(X, setdiff(seq_along(dims), axis), function(v) fft(v, inverse = inverse))
  }
  X
}

fftn <- function(x, inverse = FALSE) {
  dims <- dim(x)
  X <- x
  for (axis in seq_along(dims)) {
    n <- dims[axis]
    # bring axis to the last position
    perm <- c((1:length(dims))[-axis], axis)
    Xp <- aperm(X, perm)
    
    # flatten into matrix: rows = everything else, cols = axis
    mat <- matrix(Xp, ncol = n)
    
    # apply FFT row-wise
    mat_fft <- t(apply(mat, 1, function(row) fft(row, inverse = inverse)))
    
    # reshape back to permuted dims
    Xp <- array(mat_fft, dim = dims[perm])
    
    # restore original axis order
    invperm <- match(seq_along(dims), perm)
    X <- aperm(Xp, invperm)
  }
  X
}

ifftn <- function(x) {
  fftn(x, inverse = TRUE) / length(x)
}

# --- Smallest 3D test case ---
# A 2x2x2 cube with simple integer values
arr <- array(1:(3^4), dim = rep(3, 4))
cat("Original array:\n")
print(arr)

# Forward 3D FFT
Fft <- fftn(arr)
cat("\nFFT result (complex values):\n")
print(Fft)

# Inverse FFT
arr_rec <- ifftn(Fft)
cat("\nReconstructed array (after ifftn):\n")
print(Re(arr_rec))   # take real part (imag ~ 1e-15 noise)
