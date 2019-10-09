library(data.table)
library(ggplot)
library(ggtheme)
library(Ecdat)
library(faraway)
library(fGarch)
library(sn)

theme_set(theme_light())

setwd("D:/Projects/MSDS-RiskAnalytics/datasets")

# Theme Overrides
theme_update(plot.title = element_text(hjust = 0.5),
             axis.text.x = element_text(size = 10),
             axis.text.y = element_text(size = 10),
             axis.title = element_text(face = "bold", size = 12, colour = "steelblue4"),
             legend.position = "top", legend.title = element_blank())

#### Maximum Likelihood Lab

data(Capm, package = "Ecdat")
x <- diff(Capm$rf)
fitdistr(x, "t")

# classical t-distribution

n <- length(x)
start <- c(mean(x), sd(x), 5)
loglik_t <- function(beta) sum( - dt((x - beta[1]) / beta[2],
                                     beta[3], log = T) + log(beta[2]))
fit_t <- optim(start, loglik_t, hessian = T,
               method = "L-BFGS-B", lower = c(-1, 0.001, 1))
AIC_t <- 2 * fit_t$value + 2 * 3
BIC_t <- 2 * fit_t$value + log(n) * 3
sd_t <- sqrt(diag(solve(fit_t$hessian)))

fit_t$par
sd_t
AIC_t
BIC_t

# standardized t-distribution

loglik_std <- function(beta) sum( - dstd(x, mean = beta[1],
                                     sd = beta[2], nu = beta[3], log = T))
fit_std <- optim(start, loglik_std, hessian = T,
               method = "L-BFGS-B", lower = c(-0.1, 0.01, 2.1))
AIC_std <- 2 * fit_std$value + 2 * 3
BIC_std <- 2 * fit_std$value + log(n) * 3
sd_std <- sqrt(diag(solve(fit_std$hessian)))

fit_std$par
sd_std
AIC_std
BIC_std

# F-S skewed t-distribution

loglik_sstd <- function(beta) sum( - dsstd(x, mean = beta[1],
                                         sd = beta[2], nu = beta[3], xi = beta[4], log = T))
start <- c(mean(x), sd(x), 5, 1)
fit_sstd <- optim(start, loglik_sstd, hessian = T,
                 method = "L-BFGS-B", lower = c(-0.1, 0.01, 2.1, -2))
AIC_sstd <- 2 * fit_sstd$value + 2 * 3
BIC_sstd <- 2 * fit_sstd$value + log(n) * 3
sd_sstd <- sqrt(diag(solve(fit_sstd$hessian)))

fit_sstd$par
sd_sstd
AIC_sstd
BIC_sstd

dat <- read.csv("FlowData.csv")
dat <- dat/10000

par(mfrow = c(3,2))
x <- dat$Flow1
x1 <- sort(x)
fit1 <- sn.mple(y = x1, x = as.matrix(rep(1, length(x1))))
est1 <- cp2dp(fit1$cp, family = "SN")

plot(x1, dsn(x1, dp = est1),
     type = "l", lwd = 2, xlab = "flow",
     ylab = "flow 1 density")
d <- density(x1)
lines(d$x, d$y, lty = 2, lwd = 2)
legend(40, 0.034, c("t-model", "KDE"), lty = c(1, 2),
       lwd = c(2, 2))
n <- length(x1)
u = (1:n) / (n + 1)

plot(x1, qsn(u, dp = est1), xlab = "data",
     ylab = "skew-t quantiles", main = "Flow 1")
lmfit <- lm(qsn(c(0.25, 0.75), dp = est1) ~ quantile(x1, 
        c(0.25, 0.75)))
abline(lmfit)
