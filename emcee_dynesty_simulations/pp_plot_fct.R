


s1 = rnorm(1000, mean = 0, sd = 1)
s2 = rnorm(10000, mean = 0, sd = 1)
#s1 = rgamma(1000, shape = 7.5, scale = 1)
#s2 = rgamma(1000, shape = 7.5, scale = 1)



pp_plot_fct = function(sample1, sample2, pdf_save = FALSE, file_name = "",
                       xlab = "Quantiles sample1", ylab = "Quantiles sample2",
                       main = "P-P plot", lwd = 1){

  quantile_grid = seq(0.001, 0.999, 0.001)

  ecdf_sample1  = ecdf(sample1)
  ecdf_sample2  = ecdf(sample2)
  quant_sample1 = quantile(ecdf_sample1,  probs = quantile_grid, names = F)
  
  if (pdf_save) {
    pdf(paste(file_name, "pp_plot_tmp.pdf")) 
    plot(quantile_grid, ecdf_sample2(quant_sample1), 
         xlim = c(0,1), ylim = c(0,1),
         xlab = xlab, ylab = ylab,
         main = main)
    lines(c(0,1), c(0,1), col = "red", lwd = lwd)
    dev.off()
  } else {
    plot(quantile_grid, ecdf_sample2(quant_sample1), 
         xlim = c(0,1), ylim = c(0,1),
         xlab = xlab, ylab = ylab,
         main = main)
    lines(c(0,1), c(0,1), col = "red", lwd = lwd)
  }
}

# pdf is stored in your current qorking directory
pp_plot_fct(s1, s2)


  

