myplot = function(x, yN, seN, yS, seS, xlabel, ylabel, xl, yl) {
  


  ylabN = c("Proportion of species with",paste0(ylabel," abundance"))
  ylabS = c("Proportion of species with",paste0(ylabel," richness"))

  
  plot(x, yN,
       xlab = xlabel, ylab = ylabN,
       xlim = xl, ylim = yl,
       pch = 19)
  
  
 
  zero.se <- sum(seN == 0 | is.na(seN))
  
  if (zero.se == 0) {
    lm1 <- lm(yN ~ x, weights = 1 / seN^2)
    arrows(x0 = x, y0 = yN - seN, x1 = x, y1 = yN + seN,
           angle = 90, code = 3, length = 0.04, col = "black", lwd = 1.4)
  } else {
    lm1 <- lm(yN ~ x)
  }
  
  abline(lm1)
  
  h = summary(lm1)
  p = h$coefficients[2, 4]
  r2 = h$r.squared
  
  if (p < 1e-4) {
    title(paste0('c) R² = ', round(r2, 2), ', p < 1e-4'), cex.main = 1)
  } else {
    title(paste0('c) R² = ', round(r2, 2), ', p = ', round(p, 4)), cex.main = 1)
  }
  
  ### richness
  plot(x, yS,
       xlab = xlabel, ylab = ylabS,
       xlim = xl, ylim = yl,
       pch = 19)
  
  
  zero.se <- sum(seS == 0 | is.na(seS))
  
  if (zero.se == 0) {
    lm1 <- lm(yS ~ x, weights = 1 / seS^2)
    arrows(x0 = x, y0 = yS - seS, x1 = x, y1 = yS + seS,
           angle = 90, code = 3, length = 0.04, col = "black", lwd = 1.4)
  } else {
    lm1 <- lm(yS ~ x)
  }
  
  abline(lm1)
  
  h = summary(lm1)
  p = h$coefficients[2, 4]
  r2 = h$r.squared
  
  if (p < 1e-4) {
    title(paste0('d) R² = ', round(r2, 2), ', p < 1e-4'), cex.main = 1)
  } else {
    title(paste0('d) R² = ', round(r2, 2), ', p = ', round(p, 4)), cex.main = 1)
  }
  
}
