library(data.table)
library(xts)
library(ggplot2)
library(ggthemes)

theme_set(theme_light())

theme_set(theme_light())

# Theme Overrides
theme_update(axis.text.x = element_text(size = 10),
             axis.text.y = element_text(size = 10),
             plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "darkgreen"),
             axis.title = element_text(face = "bold", size = 12, colour = "steelblue4"),
             plot.subtitle = element_text(face = "bold", size = 8, colour = "darkred"),
             legend.title = element_text(size = 12, color = "darkred", face = "bold"),
             legend.position = "right", legend.title.align=0.5,
             panel.border = element_rect(linetype = "solid", 
                                         colour = "lightgray"), 
             plot.margin = unit(c( 0.1, 0.1, 0.1, 0.1), "inches"))

data.dir <- "D:/Projects/MSDS-RiskAnalytics/datasets/"

setwd("D:/Projects/MSDS-RiskAnalytics/Module_08")

################################################################
########## Code for figure 16.1   ##############################
################################################################

mu1 = 0.14
mu2 = 0.08
sig1 = 0.2
sig2 = 0.15
rho = 0
rf = 0.06
w = seq(0, 1, len = 500)
means = 0.08 + 0.06 * w
var = sig1^2 * w^2 + sig2^2 * (1 - w)^2
risk = sqrt(var)
ind = !(risk > min(risk))
ind2 = (means > means[ind])
wt = 0.693
meant = 0.08 + 0.06 * wt
riskt = sqrt(sig1^2 * wt^2 + sig2^2 * (1 - wt)^2)

wp = 0.475
meanp = 0.08 + 0.06 * wp
riskp = sqrt(sig1^2 * wp^2 + sig2^2 * (1 - wp)^2)


pdf("portfolioNew.pdf", width = 6, height = 5)
plot(risk, means, type = "l", lwd = 1, xlim = c(0, 0.21), ylim = c(0.0575, 0.145))
abline(v = 0)
lines(risk[ind2], means[ind2], type = "l", lwd = 5, col = "red")
lines( c(0, riskt), c(rf, meant), col = "blue", lwd =2)
lines(c(0,riskp), c(rf,meanp), col = "purple", lwd = 2, lty = 2)
text(riskt, meant, "T", cex = 1.2)
text(sig1, mu1, "R1", cex = 1.2)
text(sig2, mu2, "R2", cex = 1.2)
text(0, rf, "F", cex = 1.2)
text(riskp, meanp, "P", cex = 1.2)
text(min(risk), means[ind], "MV", cex = 1.2)
graphics.off()

