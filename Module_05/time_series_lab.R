library(Ecdat)

data(Mishkin, package = "Ecdat")

y <- as.vector(Mishkin[, 1])

par(mfrow = c(1,2))
acf(y)
acf(diff(y))

