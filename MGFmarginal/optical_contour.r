library(LaplacesDemon)
library(ggplot2)
library(viridis)
library(plotly)  # For 3D visualization

tic()
source("mgf_mle.R")
toc()

############################
# contour on (lalpha, lbeta) and (lmu, ltheta)
############################
lalpha_eval <- seq(olgammamle$par[1] - 10, olgammamle$par[1] + 1.5, length.out = 100)
lbeta_eval <- seq(olgammamle$par[2] - 12, olgammamle$par[2] + 2.5, length.out = 100)

lmu_eval <- seq(olgamma_mean - 2, olgamma_mean + 2, length.out = 100)
ltheta_eval <- seq(olgamma_var - 2, olgamma_var + 2, length.out = 100)


lalpha_loglik_fun <- function(lbeta) {
  sapply(lalpha_eval, FUN = function(lalpha) loglik_lme(c(lalpha, lbeta), oSData, bkgdalpha = bkgdShape, bkgdbeta = bkgdRates, Time = oData$Segments$Time))
}
lalphabeta_loglik <- sapply(lbeta_eval,
                            FUN = lalpha_loglik_fun)
lmu_loglik_fun <- function(ltheta) {
  sapply(lmu_eval, FUN = function(lmu) loglik_lme_meanvar(c(lmu, ltheta), oSData, bkgdalpha = bkgdShape, bkgdbeta = bkgdRates, Time = oData$Segments$Time))
}
lmutheta_loglik <- sapply(ltheta_eval,
                          FUN = lmu_loglik_fun)

contour(lmu_eval, ltheta_eval, exp(lmutheta_loglik - max(lmutheta_loglik)), nlevels = 10)
points(olgamma_mean, olgamma_var, col ="red", pch = 19)


# Create a data frame from your matrices
contour_data <- expand.grid(lalpha = lalpha_eval, lbeta = lbeta_eval)
contour_data$loglik <- as.vector(lalphabeta_loglik)
saveRDS(contour_data, "optical_shape_rate_contour.RDS")

pdf(file = "optical_shape-rate_contour.pdf", width = 5, height = 5)
# Create the ggplot
ggbins <- 100
ggplot(contour_data, aes(x = lalpha, y = lbeta, z = loglik)) +
  geom_contour_filled(aes(fill = after_stat(level_mid)), bins = ggbins, alpha = 0.8) +
  geom_contour(bins = ggbins, color = "white", alpha = 0.1, linewidth = 0.3) +
  #xlim(-3.5, -2.5) +
  #ylim(22, 23) +
  # MLE point as red cross
  geom_point(aes(x = olgammamle$par[1], y = olgammamle$par[2]),
             color = "red", size = 3, shape = 4) +
  # MLE text label
  geom_text(aes(x = olgammamle$par[1], y = olgammamle$par[2], 
                label = "MLE", color = "red"), vjust = 1.5, size = 4.5) +
  scale_fill_viridis_c(option = "plasma", name = "log likelihood") +
  scale_color_identity() +  # Ensures text uses exact color specified
  labs(title = "Likelihood Contours",
       subtitle = "isolated & binary-overlapping optical sources",
       x = "log shape",
       y = "log rate") +
  theme_minimal() +
  theme(legend.position = "right")

ggplot(contour_data, aes(x = lalpha, y = lbeta, z = loglik)) +
  geom_contour_filled(aes(fill = after_stat(level_mid)), bins = ggbins, alpha = 0.8) +
  geom_contour(bins = ggbins, color = "white", alpha = 0.1, linewidth = 0.3) +
  # MLE point as red cross
  geom_point(aes(x = olgammamle$par[1], y = olgammamle$par[2]),
             color = "red", size = 3, shape = 4) +
  # MLE text label
  geom_text(aes(x = olgammamle$par[1], y = olgammamle$par[2], 
                label = "MLE", color = "red"), vjust = 1.5, size = 4.5) +
  scale_fill_viridis_c(option = "viridis", name = "log likelihood") +
  scale_color_identity() +  # Ensures text uses exact color specified
  labs(title = "Likelihood Contours",
       subtitle = "isolated & binary-overlapping optical sources",
       x = "log shape",
       y = "log rate") +
  theme_minimal() +
  theme(legend.position = "right")

plot(lalpha_eval, rowMeans(lalphabeta_loglik), type = "l", 
     main = "Marginal likelihood", sub = "with uniform weight", xlab = "log(alpha)", ylab = "log likelihood")
plot(lbeta_eval, colMeans(lalphabeta_loglik), type = "l", 
     main = "Marginal likelihood", sub = "with uniform weight", xlab = "log(beta)", ylab = "log likelihood")
dev.off()

###############
# 3D plot on (lmu, ltheta)
################
# Create a data frame from your matrices
contour_meanvar <- expand.grid(lmu = lmu_eval, ltheta = ltheta_eval)
contour_meanvar$loglik <- as.vector(lmutheta_loglik)
saveRDS(contour_meanvar, "opticalmu_theta_contour.RDS")
contour_meanvar <- readRDS("optical_mu_theta_contour.RDS")

########
# priors
########
# cauchy
lcauchy_mu <- dcauchy(contour_meanvar$lmu, location = olgamma_mean, scale = exp(olgamma_mean) * 100, log = TRUE)
lcauchy_theta <- dcauchy(contour_meanvar$ltheta, location = olgamma_var, scale = exp(olgamma_var) * 100, log = TRUE)
contour_meanvar$cauchy_post <- lcauchy_mu + lcauchy_theta + contour_meanvar$loglik

# t4
lt4_mu <- dt(contour_meanvar$lmu, df = 4, log = TRUE)
lt4_theta <- dt(contour_meanvar$ltheta, df = 4, log = TRUE)
contour_meanvar$t4_post <- lt4_mu + lt4_theta + contour_meanvar$loglik

# t1
lt1_mu <- dt(contour_meanvar$lmu, df = 1, log = TRUE)
lt1_theta <- dt(contour_meanvar$ltheta, df = 1, log = TRUE)
contour_meanvar$t1_post <- lt1_mu + lt1_theta + contour_meanvar$loglik

# normal
lnorm_mu <- dnorm(contour_meanvar$lmu, mean = olgamma_mean, sd = sqrt(exp(olgamma_mean) * 100), log = TRUE)
lnorm_theta <- dnorm(contour_meanvar$ltheta, mean = olgamma_var, sd = sqrt(exp(olgamma_var) * 100), log = TRUE)
contour_meanvar$norm_post <- lnorm_mu + lnorm_theta + contour_meanvar$loglik


pdf(file = "optical-mean-variance_contour.pdf", width = 5, height = 5)
# Create the ggplot
ggbins <- 25
ggplot(contour_meanvar, aes(x = lmu, y = ltheta, z = exp(loglik - max(loglik)))) +
  geom_contour_filled(aes(fill = after_stat(level_mid)), bins = ggbins, alpha = 0.8) +
  #geom_contour(bins = ggbins, color = "white", alpha = 1e-6, linewidth = 0.3) +
  xlim(-25.75, -25) +
  ylim(-49, -47) +
  # MLE point as red cross
  geom_point(aes(x = olgamma_mean, y = olgamma_var),
             color = "red", size = 3, shape = 4) +
  # MLE text label
  geom_text(aes(x = olgamma_mean, y = olgamma_var, 
                label = "MLE", color = "red"), vjust = 1.5, size = 4.5) +
  scale_fill_viridis_c(option = "plasma", name = "scaled likelihood") +
  scale_color_identity() +  # Ensures text uses exact color specified
  labs(title = "Likelihood Contours",
       subtitle = "isolated & binary-overlapping optical sources",
       x = "log mean",
       y = "log variance") +
  theme_minimal() +
  theme(legend.position = "right")

ggplot(contour_meanvar, aes(x = lmu, y = ltheta, z = exp(loglik - max(loglik)))) +
  geom_contour_filled(aes(fill = after_stat(level_mid)), bins = ggbins, alpha = 0.8) +
  #geom_contour(bins = ggbins, color = "white", alpha = 0.1, linewidth = 0.3) +
  xlim(-25.75, -25) +
  ylim(-49, -47) +
  # MLE point as red cross
  geom_point(aes(x = olgamma_mean, y = olgamma_var),
             color = "red", size = 3, shape = 4) +
  # MLE text label
  geom_text(aes(x = olgamma_mean, y = olgamma_var, 
                label = "MLE", color = "red"), vjust = 1.5, size = 4.5) +
  scale_fill_viridis_c(option = "viridis", name = "scaled likelihood") +
  scale_color_identity() +  # Ensures text uses exact color specified
  labs(title = "Likelihood Contours",
       subtitle = "isolated & binary-overlapping optical sources",
       x = "log mean",
       y = "log variance") +
  theme_minimal() +
  theme(legend.position = "right")


dev.off()



# other contours
ggplot(contour_meanvar, aes(x = lmu, y = ltheta, z = cauchy_post)) +
  geom_contour_filled(aes(fill = after_stat(level_mid)), bins = ggbins, alpha = 0.8) +
  geom_contour(bins = ggbins, color = "white", alpha = 0.1, linewidth = 0.3) +
  # MLE point as red cross
  geom_point(aes(x = olgamma_mean, y = olgamma_var),
             color = "red", size = 3, shape = 4) +
  # MLE text label
  geom_text(aes(x = olgamma_mean, y = olgamma_var, 
                label = "MLE", color = "red"), vjust = 1.5, size = 4.5) +
  scale_fill_viridis_c(option = "viridis", name = "posterior kernel") +
  scale_color_identity() +  # Ensures text uses exact color specified
  labs(title = "Cauchy(mle, 100*mle) prior",
       x = "log mean",
       y = "log variance") +
  theme_minimal() +
  theme(legend.position = "right")

ggplot(contour_meanvar, aes(x = lmu, y = ltheta, z = t4_post)) +
  geom_contour_filled(aes(fill = after_stat(level_mid)), bins = ggbins, alpha = 0.8) +
  geom_contour(bins = ggbins, color = "white", alpha = 0.1, linewidth = 0.3) +
  # MLE point as red cross
  geom_point(aes(x = olgamma_mean, y = olgamma_var),
             color = "red", size = 3, shape = 4) +
  # MLE text label
  geom_text(aes(x = olgamma_mean, y = olgamma_var, 
                label = "MLE", color = "red"), vjust = 1.5, size = 4.5) +
  scale_fill_viridis_c(option = "viridis", name = "posterior kernel") +
  scale_color_identity() +  # Ensures text uses exact color specified
  labs(title = "t(df=4) prior",
       x = "log mean",
       y = "log variance") +
  theme_minimal() +
  theme(legend.position = "right")

ggplot(contour_meanvar, aes(x = lmu, y = ltheta, z = t1_post)) +
  geom_contour_filled(aes(fill = after_stat(level_mid)), bins = ggbins, alpha = 0.8) +
  geom_contour(bins = ggbins, color = "white", alpha = 0.1, linewidth = 0.3) +
  # MLE point as red cross
  geom_point(aes(x = olgamma_mean, y = olgamma_var),
             color = "red", size = 3, shape = 4) +
  # MLE text label
  geom_text(aes(x = olgamma_mean, y = olgamma_var, 
                label = "MLE", color = "red"), vjust = 1.5, size = 4.5) +
  scale_fill_viridis_c(option = "viridis", name = "posterior kernel") +
  scale_color_identity() +  # Ensures text uses exact color specified
  labs(title = "t(df=1) prior",
       x = "log mean",
       y = "log variance") +
  theme_minimal() +
  theme(legend.position = "right")

ggplot(contour_meanvar, aes(x = lmu, y = ltheta, z = norm_post)) +
  geom_contour_filled(aes(fill = after_stat(level_mid)), bins = ggbins, alpha = 0.8) +
  geom_contour(bins = ggbins, color = "white", alpha = 0.1, linewidth = 0.3) +
  # MLE point as red cross
  geom_point(aes(x = olgamma_mean, y = olgamma_var),
             color = "red", size = 3, shape = 4) +
  # MLE text label
  geom_text(aes(x = olgamma_mean, y = olgamma_var, 
                label = "MLE", color = "red"), vjust = 1.5, size = 4.5) +
  scale_fill_viridis_c(option = "viridis", name = "posterior kernel") +
  scale_color_identity() +  # Ensures text uses exact color specified
  labs(title = "normal(mle, 100*mle) prior",
       x = "log mean",
       y = "log variance") +
  theme_minimal() +
  theme(legend.position = "right")


plot(lmu_eval, rowMeans(lmutheta_loglik), type = "l", 
     main = "Marginal likelihood", sub = "with uniform weight", xlab = "log(mu)", ylab = "log likelihood")
plot(ltheta_eval, colMeans(lmutheta_loglik), type = "l", 
     main = "Marginal likelihood", sub = "with uniform weight", xlab = "log(theta)", ylab = "log likelihood")


# interactive 3D plot for (lmu, ltheta)
# Convert your data to a matrix format (required for surface plots)
z_matrix <- matrix(contour_meanvar$loglik, 
                   nrow = length(unique(contour_meanvar$lmu)),
                   ncol = length(unique(contour_meanvar$ltheta)))
z_matrix <- lmutheta_loglik[nrow(lmutheta_loglik):1, ncol(lmutheta_loglik):1]
# Create the 3D surface plot
plot_ly() %>%
  add_surface(
    x = unique(contour_meanvar$lmu),
    y = unique(contour_meanvar$ltheta),
    z = z_matrix,
    colorscale = list(c(0, 0.5, 1), c("#0D0887", "#CC4678", "#F0F921")),  # Plasma palette
    opacity = 0.8,
    contours = list(
      z = list(show = TRUE, color = "white", width = 0.3)  # White contour lines
    )
  ) %>%
  add_markers(
    x = olgamma_mean,
    y = olgamma_var,
    z = olgammamle$value,
    marker = list(color = "red", size = 5, symbol = "x"),
    name = "MLE"
  ) %>%
  add_annotations(
    x = olgamma_mean,
    y = olgamma_var,
    z = olgammamle$value,
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
      SData,
      lalpha = lalpha,
      lbeta = lbeta,
      Time = xData$Segments$Time,
      bkgdalpha = bkgdShape,
      bkgdbeta = bkgdRates
    )
    return(l)
  }
}
logker_cube_meanvar(c(0.6, 0.6), oSData)
lmu_cube <- seq(olgamma_mean - 15, olgamma_mean + 2, length.out = 100)
ltheta_eval <- seq(olgamma_var - 5, olgamma_var + 45, length.out = 100)


plot(seq(-1, 1, by = 1e-2), 
     dcauchy(seq(-1, 1, by = 1e-2), location = 1.2e-6, scale = 1.2e-6 * 100), type = "l",
     xlab = "parameter", ylab = "density", main = "Cauchy(1.2e-6, 1.2e-4)")
plot(seq(-100, 100, by = 1), 
     dcauchy(seq(-100, 100, by = 1), location = olgamma_mean, scale = exp(olgamma_mean) * 100), type = "l",
     xlab = "parameter", ylab = "density", main = "Cauchy(0, 1e6)")
plot(seq(-100, 100, by = 1), 
     dnorm(seq(-100, 100, by = 1), mean = log(1.2e-6), sd = 1), type = "l",
     xlab = "parameter", ylab = "density", main = "Cauchy(0, 1e6)")
