library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(GGally)
library(ggthemes)
library(scales)
library(reshape2)
library(skimr)
library(gridExtra)

#####################################################################
######################### Computation 6 #############################
#####################################################################

theme_set(theme_sjplot())

path.data <- "D:/Projects/MSDS-RiskAnalytics/datasets"
setwd(path.data)

dat <- read.csv("Stock_bond.csv", header = T)

names(dat)
attach(dat)

par(mfrow = c(1, 2))

plot(GM_AC, type = "l")
plot(F_AC, type = "l")

n <- dim(dat)[1]
GMReturn <- GM_AC[-1] / GM_AC[-n] - 1
FReturn <- F_AC[-1] / F_AC[-n] - 1

par(mfrow = c(1, 1))
plot(GMReturn, FReturn)

# Problems

# 1.)
# Do GM and Ford returns seem positively correlated? Do you notice any outliying returns? 

cor(GMReturn, FReturn) # .61

# If "yes," do outlying GM returns seem to occur with outliying Ford returns?

logGMReturn <- diff(log(GM_AC))

plot(GMReturn, logGMReturn)
cor(GMReturn, logGMReturn)

# MSFT / MRK

MSReturn <- MSFT_AC[-1] / MSFT_AC[-n] - 1
logMSReturn <- diff(log(MSFT_AC))

MRKReturn <- MRK_AC[-1] / MRK_AC[-n] - 1
logMRKReturn <- diff(log(MRK_AC))

plot(MSReturn, MRKReturn)

cor(MSReturn, MRKReturn)

plot(MSReturn, logMSReturn)

# Simulations

# Suppose a hedge fund owns $1,000,000 of stock and used $50,000 of its own capital and $950,000 in borrowed money for the purchase.
# Suppose that if the value of the stock falls below $950,000 at the end of any trading day, then the hedge fund will sell all the stock and repay the loan.
# This will wipe out its $50,000 investment.
# The hedge fund is said to be leveraged 20:1 since its position is 220 times the amount of its own capital invested.

# Suppose that the daily log returns on the stock have a mean of 0.05/year and a standard deviation of 0.23/year. 
# These can be converted to rates per trading day by dividing by 253 and sqrt(253), respectively.

# What is the probability that the value of the stock will be below $950,000 at the close of at least one of the next 45 trading days?

niter <- 1e5 # number of iterations
set.seed(2009) # reproducible

# setup storage
outcomes <- list(below = rep(0, niter),
                  pnl = rep(0, niter),
                  ret = rep(0, niter))

seed.capital <- 5e4
initial.investment <- log(1e6)
profit.threshold <- log(1.1e6)
loss.threshold <- log(9.5e5)

simulate_market <- function(days) {
  # generate random returns for N days
  r <- rnorm(days, mean = 0.05 / 253,
            sd = 0.23 / sqrt(253))

  cumsum(r) # return the final log price after N days.
}

for (i in 1:niter) {

  logPrice = initial.investment + simulate_market(45) # simulate 45 trading days.

  minlogP = min(logPrice) # miniumum price over next 45 days

  outcomes$below[i] = as.numeric(minlogP < loss.threshold)
}

print(paste0("Probability the value of the stock is below $950,000 at least one of next 45 sessions: ", mean(outcomes$below) * 100, "%"))

# reset seed for next simulation.

set.seed(2009) # reproducible

for (i in 1:niter) {

  logPrice = initial.investment + simulate_market(100) # simulate 100 trading days.

  suppressWarnings({
    # ignore Inf returned if condition not meet.
    profit.day <- min(which(logPrice >= profit.threshold))
    loss.day <- min(which(logPrice <= loss.threshold))
  })

  # What was the exit condition of the position, hince the final price of the stock?
  daysOpen <- ifelse(profit.day == Inf && loss.day == Inf, length(logPrice),
                       min(profit.day, loss.day))

  # p&l = ending value - initial investment
  outcomes$pnl[i] <- exp(logPrice[daysOpen]) - exp(initial.investment)
  outcomes$ret[i] <- outcomes$pnl[i] / seed.capital
}

sum(outcomes$pnl > 0) / length(outcomes$pnl) # Probability of Profit
sum(outcomes$pnl < 0) / length(outcomes$pnl) # Probability of Loss
mean(outcomes$pnl) # Expected P&L


# execution.type <- ifelse(profit.day == Inf && loss.day == Inf, "Market",
#  ifelse(profit.day < loss.day, "Profit", "Loss"))

# execution.type

# priceIndex <- ifelse(execution.type == "Market", length(logPrice), min(profit.day, loss.day))
