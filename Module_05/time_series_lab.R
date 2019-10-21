library(Ecdat)

data(Mishkin, package = "Ecdat")

y <- as.vector(Mishkin[, 1])

par(mfrow = c(1,2))
acf(y)
acf(diff(y))

Box.test(diff(y), lag = 10, type = "Ljung-Box")

data(bmw, package = "evir")

Box.test(bmw, lag = 5, type = "Ljung-Box")

fitAR1 <- arima(bmw, order = c(1, 0, 0))
summary(fitAR1)

acf(bmw)

Box.test(residuals(fitAR1), lag = 5, type = "Ljung-Box", fitdf = 1)

fit <- arima(y, order = c(1, 0, 0))
Box.test(fit$residuals, type = "Ljung", lag = 24, fitdf = 1)

library(forecast)

auto.arima(diff(y), max.p = 20, max.q = 0, d = 0, ic = "aic")

