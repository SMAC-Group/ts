library(astsa)
library(mgcv)
library(simts)

### Chapter 1
plot(globtemp, type="o", ylab="Global Temperature Deviations")

time <- 1:length(globtemp)
fit <- gam(globtemp ~ s(time))
resids <- resid(fit)

plot(resids, type="o", ylab="Global Temperature Deviations")

time <- 1:100
set.seed(9)
y.ind <- cumsum(rep(0.01, 100)) + rnorm(100)
fit.ind <- lm(y.ind ~ time)
summary(fit.ind)
y.dep <- cumsum(rep(0.01, 100)) + arima.sim(n = 100, list(ar = c(0.8897, -0.4858)))
fit.dep <- lm(y.dep ~ time)
summary(fit.dep)

time_jj <- 1:84
quarters <- as.factor(rep(1:4, 21))
fit_jj1 <- lm(as.vector(jj) ~ time_jj)
plot(time_jj, as.vector(jj), type="l", ylab="Quarterly Earnings per Share")
abline(fit_jj1, col = "red")

fit_jj2 <- lm(as.vector(jj) ~ quarters + time_jj + I(time_jj^2))
plot(time_jj, as.vector(jj), type="l", ylab="Quarterly Earnings per Share")
lines(time_jj, fit_jj2$fitted.values, col = "red")

fit_jj3 <- gam(as.vector(jj) ~ quarters + s(time_jj))
plot(time_jj, as.vector(jj), type="l", ylab="Quarterly Earnings per Share")
lines(time_jj, fit_jj3$fitted.values, col = "red")

plot(fit_jj3$residuals, type = "l")
