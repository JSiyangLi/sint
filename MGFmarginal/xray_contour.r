library(LaplacesDemon)
library(ggplot2)
library(viridis)
library(plotly)  # For 3D visualization
library(tictoc)

tic()
source("mgf_mle.R")
toc()

############################
# contour on (lalpha, lbeta) and (lmu, ltheta)
############################
xlalpha_eval <- seq(xlgammamle$par[1] - 10, xlgammamle$par[1] + 1.5, length.out = 100)
xlbeta_eval <- seq(xlgammamle$par[2] - 12, xlgammamle$par[2] + 2.5, length.out = 100)

xlmu_eval <- seq(xlgamma_mean - 2, xlgamma_mean + 2, length.out = 100)
xltheta_eval <- seq(xlgamma_var - 2, xlgamma_var + 2, length.out = 100)


xlalpha_loglik_fun <- function(lbeta) {
  sapply(xlalpha_eval, FUN = function(lalpha) loglik_lme(c(lalpha, lbeta), xSData, 
                                                         bkgdalpha = bkgdShape, bkgdbeta = bkgdRates, Time = xData$Segments$Time))
}
xlalphabeta_loglik <- sapply(xlbeta_eval,
                            FUN = xlalpha_loglik_fun)
xlmu_loglik_fun <- function(ltheta) {
  sapply(xlmu_eval, FUN = function(lmu) loglik_lme_meanvar(c(lmu, ltheta), xSData,
                                                           bkgdalpha = bkgdShape, bkgdbeta = bkgdRates, Time = xData$Segments$Time))
}
xlmutheta_loglik <- sapply(xltheta_eval,
                          FUN = xlmu_loglik_fun)

contour(xlmu_eval, xltheta_eval, xlmutheta_loglik, nlevels = 100)
points(xlgamma_mean, xlgamma_var, col ="red", pch = 19)


# Create a data frame from your matrices
xcontour_data <- expand.grid(lalpha = xlalpha_eval, lbeta = xlbeta_eval)
xcontour_data$loglik <- as.vector(xlalphabeta_loglik)
saveRDS(xcontour_data, "xray_shape_rate_contour.RDS")

pdf(file = "xray_shape-rate_contour.pdf", width = 5, height = 5)
# Create the ggplot
ggbins <- 100
ggplot(xcontour_data, aes(x = lalpha, y = lbeta, z = loglik)) +
  geom_contour_filled(aes(fill = after_stat(level_mid)), bins = ggbins, alpha = 0.8) +
  geom_contour(bins = ggbins, color = "white", alpha = 0.1, linewidth = 0.3) +
  # MLE point as red cross
  geom_point(aes(x = xlgammamle$par[1], y = xlgammamle$par[2]),
             color = "red", size = 3, shape = 4) +
  # MLE text label
  geom_text(aes(x = xlgammamle$par[1], y = xlgammamle$par[2], 
                label = "MLE", color = "red"), vjust = 1.5, size = 4.5) +
  scale_fill_viridis_c(option = "plasma", name = "log likelihood") +
  scale_color_identity() +  # Ensures text uses exact color specified
  labs(title = "Likelihood Contours",
       subtitle = "of all Xray sources",
       x = "log shape",
       y = "log rate") +
  theme_minimal() +
  theme(legend.position = "right")

ggplot(xcontour_data, aes(x = lalpha, y = lbeta, z = loglik)) +
  geom_contour_filled(aes(fill = after_stat(level_mid)), bins = ggbins, alpha = 0.8) +
  geom_contour(bins = ggbins, color = "white", alpha = 0.1, linewidth = 0.3) +
  # MLE point as red cross
  geom_point(aes(x = xlgammamle$par[1], y = xlgammamle$par[2]),
             color = "red", size = 3, shape = 4) +
  # MLE text label
  geom_text(aes(x = xlgammamle$par[1], y = xlgammamle$par[2], 
                label = "MLE", color = "red"), vjust = 1.5, size = 4.5) +
  scale_fill_viridis_c(option = "viridis", name = "log likelihood") +
  scale_color_identity() +  # Ensures text uses exact color specified
  labs(title = "Likelihood Contours",
       subtitle = "of all Xray sources",
       x = "log shape",
       y = "log rate") +
  theme_minimal() +
  theme(legend.position = "right")

plot(xlalpha_eval, rowMeans(xlalphabeta_loglik), type = "l", 
     main = "Marginal likelihood", sub = "with uniform weight", xlab = "log(alpha)", ylab = "log likelihood")
plot(xlbeta_eval, colMeans(xlalphabeta_loglik), type = "l", 
     main = "Marginal likelihood", sub = "with uniform weight", xlab = "log(beta)", ylab = "log likelihood")
dev.off()

###############
# 3D plot on (lmu, ltheta)
################
# Create a data frame from your matrices
xcontour_meanvar <- expand.grid(lmu = xlmu_eval, ltheta = xltheta_eval)
xcontour_meanvar$loglik <- as.vector(xlmutheta_loglik)
saveRDS(xcontour_meanvar, "xray_mu_theta_contour.RDS")
xcontour_meanvar <- readRDS("xray_mu_theta_contour.RDS")

########
# priors
########
# cauchy
lcauchy_mu <- dcauchy(xcontour_meanvar$lmu, location = xlgamma_mean, scale = exp(xlgamma_mean) * 100, log = TRUE)
lcauchy_theta <- dcauchy(xcontour_meanvar$ltheta, location = xlgamma_var, scale = exp(xlgamma_var) * 100, log = TRUE)
xcontour_meanvar$cauchy_post <- lcauchy_mu + lcauchy_theta + xcontour_meanvar$loglik

# t4
lt4_mu <- dt(xcontour_meanvar$lmu, df = 4, log = TRUE)
lt4_theta <- dt(xcontour_meanvar$ltheta, df = 4, log = TRUE)
xcontour_meanvar$t4_post <- lt4_mu + lt4_theta + xcontour_meanvar$loglik

# t1
lt1_mu <- dt(xcontour_meanvar$lmu, df = 1, log = TRUE)
lt1_theta <- dt(xcontour_meanvar$ltheta, df = 1, log = TRUE)
xcontour_meanvar$t1_post <- lt1_mu + lt1_theta + xcontour_meanvar$loglik

# normal
lnorm_mu <- dnorm(xcontour_meanvar$lmu, mean = xlgamma_mean, sd = sqrt(exp(xlgamma_mean) * 100), log = TRUE)
lnorm_theta <- dnorm(xcontour_meanvar$ltheta, mean = xlgamma_var, sd = sqrt(exp(xlgamma_var) * 100), log = TRUE)
xcontour_meanvar$norm_post <- lnorm_mu + lnorm_theta + xcontour_meanvar$loglik


pdf(file = "xray_mean-variance_contour.pdf", width = 5, height = 5)
# Create the ggplot
ggbins <- 100
ggplot(xcontour_meanvar, aes(x = lmu, y = ltheta, z = exp(loglik - max(loglik)))) +
  geom_contour_filled(aes(fill = after_stat(level_mid)), bins = ggbins, alpha = 0.8) +
  #geom_contour(bins = ggbins, color = "white", alpha = 1e-6, linewidth = 0.3) +
  xlim(-23, -22) +
  ylim(-46, -44) +
  # MLE point as red cross
  geom_point(aes(x = xlgamma_mean, y = xlgamma_var),
             color = "red", size = 3, shape = 4) +
  # MLE text label
  geom_text(aes(x = xlgamma_mean, y = xlgamma_var, 
                label = "MLE", color = "red"), vjust = 1.5, size = 4.5) +
  scale_fill_viridis_c(option = "plasma", name = "scaled likelihood") +
  scale_color_identity() +  # Ensures text uses exact color specified
  labs(title = "Likelihood Contours",
       subtitle = "of all Xray sources",
       x = "log mean",
       y = "log variance") +
  theme_minimal() +
  theme(legend.position = "right")

ggplot(xcontour_meanvar, aes(x = lmu, y = ltheta, z = exp(loglik - max(loglik)))) +
  geom_contour_filled(aes(fill = after_stat(level_mid)), bins = ggbins, alpha = 0.8) +
  #geom_contour(bins = ggbins, color = "white", alpha = 0.1, linewidth = 0.3) +
  xlim(-23, -22) +
  ylim(-46, -44) +
  # MLE point as red cross
  geom_point(aes(x = xlgamma_mean, y = xlgamma_var),
             color = "red", size = 3, shape = 4) +
  # MLE text label
  geom_text(aes(x = xlgamma_mean, y = xlgamma_var, 
                label = "MLE", color = "red"), vjust = 1.5, size = 4.5) +
  scale_fill_viridis_c(option = "viridis", name = "scaled likelihood") +
  scale_color_identity() +  # Ensures text uses exact color specified
  labs(title = "Likelihood Contours",
       subtitle = "of all Xray sources",
       x = "log mean",
       y = "log variance") +
  theme_minimal() +
  theme(legend.position = "right")

plot(xlmu_eval, rowMeans(xlmutheta_loglik), type = "l", 
     main = "Marginal likelihood", sub = "with uniform weight", xlab = "log(mu)", ylab = "log likelihood")
plot(xltheta_eval, colMeans(xlmutheta_loglik), type = "l", 
     main = "Marginal likelihood", sub = "with uniform weight", xlab = "log(theta)", ylab = "log likelihood")
dev.off()

ggplot(xcontour_meanvar, aes(x = lmu, y = ltheta, z = cauchy_post)) +
  geom_contour_filled(aes(fill = after_stat(level_mid)), bins = ggbins, alpha = 0.8) +
  geom_contour(bins = ggbins, color = "white", alpha = 0.1, linewidth = 0.3) +
  # MLE point as red cross
  geom_point(aes(x = xlgamma_mean, y = xlgamma_var),
             color = "red", size = 3, shape = 4) +
  # MLE text label
  geom_text(aes(x = xlgamma_mean, y = xlgamma_var, 
                label = "MLE", color = "red"), vjust = 1.5, size = 4.5) +
  scale_fill_viridis_c(option = "viridis", name = "posterior kernel") +
  scale_color_identity() +  # Ensures text uses exact color specified
  labs(title = "Cauchy(mle, 100*mle) prior",
       x = "log mean",
       y = "log variance") +
  theme_minimal() +
  theme(legend.position = "right")

ggplot(xcontour_meanvar, aes(x = lmu, y = ltheta, z = t4_post)) +
  geom_contour_filled(aes(fill = after_stat(level_mid)), bins = ggbins, alpha = 0.8) +
  geom_contour(bins = ggbins, color = "white", alpha = 0.1, linewidth = 0.3) +
  # MLE point as red cross
  geom_point(aes(x = xlgamma_mean, y = xlgamma_var),
             color = "red", size = 3, shape = 4) +
  # MLE text label
  geom_text(aes(x = xlgamma_mean, y = xlgamma_var, 
                label = "MLE", color = "red"), vjust = 1.5, size = 4.5) +
  scale_fill_viridis_c(option = "viridis", name = "posterior kernel") +
  scale_color_identity() +  # Ensures text uses exact color specified
  labs(title = "t(df=4) prior",
       x = "log mean",
       y = "log variance") +
  theme_minimal() +
  theme(legend.position = "right")

ggplot(xcontour_meanvar, aes(x = lmu, y = ltheta, z = t1_post)) +
  geom_contour_filled(aes(fill = after_stat(level_mid)), bins = ggbins, alpha = 0.8) +
  geom_contour(bins = ggbins, color = "white", alpha = 0.1, linewidth = 0.3) +
  # MLE point as red cross
  geom_point(aes(x = xlgamma_mean, y = xlgamma_var),
             color = "red", size = 3, shape = 4) +
  # MLE text label
  geom_text(aes(x = xlgamma_mean, y = xlgamma_var, 
                label = "MLE", color = "red"), vjust = 1.5, size = 4.5) +
  scale_fill_viridis_c(option = "viridis", name = "posterior kernel") +
  scale_color_identity() +  # Ensures text uses exact color specified
  labs(title = "t(df=1) prior",
       x = "log mean",
       y = "log variance") +
  theme_minimal() +
  theme(legend.position = "right")

ggplot(xcontour_meanvar, aes(x = lmu, y = ltheta, z = norm_post)) +
  geom_contour_filled(aes(fill = after_stat(level_mid)), bins = ggbins, alpha = 0.8) +
  geom_contour(bins = ggbins, color = "white", alpha = 0.1, linewidth = 0.3) +
  # MLE point as red cross
  geom_point(aes(x = xlgamma_mean, y = xlgamma_var),
             color = "red", size = 3, shape = 4) +
  # MLE text label
  geom_text(aes(x = xlgamma_mean, y = xlgamma_var, 
                label = "MLE", color = "red"), vjust = 1.5, size = 4.5) +
  scale_fill_viridis_c(option = "viridis", name = "posterior kernel") +
  scale_color_identity() +  # Ensures text uses exact color specified
  labs(title = "normal(mle, 100*mle) prior",
       x = "log mean",
       y = "log variance") +
  theme_minimal() +
  theme(legend.position = "right")



# interactive 3D plot for (lmu, ltheta)
# Convert your data to a matrix format (required for surface plots)
z_matrix <- matrix(xcontour_meanvar$loglik, 
                   nrow = length(unique(xcontour_meanvar$lmu)),
                   ncol = length(unique(xcontour_meanvar$ltheta)))
z_matrix <- xlmutheta_loglik[nrow(xlmutheta_loglik):1, ncol(xlmutheta_loglik):1]
# Create the 3D surface plot
plot_ly() %>%
  add_surface(
    x = unique(xcontour_meanvar$lmu),
    y = unique(xcontour_meanvar$ltheta),
    z = z_matrix,
    colorscale = list(c(0, 0.5, 1), c("#0D0887", "#CC4678", "#F0F921")),  # Plasma palette
    opacity = 0.8,
    contours = list(
      z = list(show = TRUE, color = "white", width = 0.3)  # White contour lines
    )
  ) %>%
  add_markers(
    x = xlgamma_mean,
    y = xlgamma_var,
    z = xlgammamle$value,
    marker = list(color = "red", size = 5, symbol = "x"),
    name = "MLE"
  ) %>%
  add_annotations(
    x = xlgamma_mean,
    y = xlgamma_var,
    z = xlgammamle$value,
    text = "MLE",
    showarrow = FALSE,
    font = list(color = "red", size = 12),
    yshift = 20
  ) %>%
  layout(
    title = "3D Log-Likelihood Surface",
    scene = list(
      xaxis = list(title = "log mean"),
      yaxis = list(title = "log variance"),
      zaxis = list(title = "log likelihood"),
      camera = list(eye = list(x = -1.5, y = -1.5, z = 0.5))  # Adjust view angle
    )
  )


#####################
# contour on cube of (lmu, ltheta)
#####################
logker_cube_meanvar <- function(lparams, SData) {
  lmu <- qcauchy(lparams[1], location = 0, scale = 100)
  ltheta <- qcauchy(lparams[2], location = 0, scale = 100)
  
  lalpha <- 2 * lmu - ltheta
  lbeta <- lmu - ltheta
  if (
    any(c(lalpha, lbeta) <= -.Machine$double.xmax) ||
    any(c(lalpha, lbeta) >= .Machine$double.xmax)
  ) {
    return(as.numeric(-1e100))
  } else {
    bkgdShape <- 1e-06 + xData$Background$cnt
    bkgdRates <- xData$Segments$Time * xData$Background$area * (1 + 1e-12)
    l <- lmgfMarg(
      xSData,
      lalpha = lalpha,
      lbeta = lbeta,
      Time = xData$Segments$Time,
      bkgdalpha = bkgdShape,
      bkgdbeta = bkgdRates
    )
    return(l)
  }
}
logker_cube_meanvar(c(0.6, 0.6), xSData)
lmu_cube <- seq(xlgamma_mean - 15, xlgamma_mean + 2, length.out = 100)
xltheta_eval <- seq(xlgamma_var - 5, xlgamma_var + 45, length.out = 100)


plot(seq(-1, 1, by = 1e-2), 
     dcauchy(seq(-1, 1, by = 1e-2), location = 1.2e-6, scale = 1.2e-6 * 100), type = "l",
     xlab = "parameter", ylab = "density", main = "Cauchy(1.2e-6, 1.2e-4)")
plot(seq(-100, 100, by = 1), 
     dcauchy(seq(-100, 100, by = 1), location = xlgamma_mean, scale = exp(xlgamma_mean) * 100), type = "l",
     xlab = "parameter", ylab = "density", main = "Cauchy(0, 1e6)")
plot(seq(-100, 100, by = 1), 
     dnorm(seq(-100, 100, by = 1), mean = log(1.2e-6), sd = 1), type = "l",
     xlab = "parameter", ylab = "density", main = "Cauchy(0, 1e6)")