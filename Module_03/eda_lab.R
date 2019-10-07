library(data.table)
library(ggplot)
library(ggtheme)
library(Ecdat)

theme_set(theme_light())

# Theme Overrides
theme_update(plot.title = element_text(hjust = 0.5),
             axis.text.x = element_text(size = 10),
             axis.text.y = element_text(size = 10),
             axis.title = element_text(face = "bold", size = 12, colour = "steelblue4"),
             legend.position = "top", legend.title = element_blank())

data(SP500, package = "Ecdat")

SPreturn <- SP500$r500
n <- length(SPreturn)
year_SP <- 1981 + (1:n) * (1991.25 - 1981) / n

# Visualize S&P500 log returns
plot(year_SP, SPreturn, main = "S&P 500 daily log returns",
     xlab = "year", type = "l", ylab = "log return")

ggplot(data.table(date = year_SP, return = SPreturn, direction = ifelse(SPreturn >= 0, "red", "green")), 
       aes(date, return, color = direction)) +
  geom_bar(stat = 'identity') +
  theme(legend.position = "none")


set.seed("991155")

edf_norm <- ecdf(rnorm(150))
pdf("normalcdfplot.pdf", width = 6, height = 5)
par(mfrow = c(1, 1))
plot(edf_norm, verticals = T, do.p = F, main = "EDF and CDF")
tt = seq(from = -3, to = 3, by = 0.01)
lines(tt, pnorm(tt),lty = 2, lwd = 2, col = "red")))
legend(1.5, 0.2, c("EDF", "CDF"), lty = c(1, 2),
       lwd = c(1.5, 2), col = c("black", "red"))
graphics.off()