# p.1Dis <- as.matrix(read.csv2(file.path("PSA", "Output", "PSA", "P_1Dint_scale.csv")))[, -1]
# dimnames(p.1Dis) <- list(c("Intercept", "Log_scale"))

# paramDSA <- function(x, min, max, steps = 7){
#   xout <- seq(min, max, length.out = steps*2 + 1)
#   xout <- xout[-(steps + 1)]
#   xout <- c(x, xout)
#   return(xout)
# }

# Range of intercepts and log_scales: 95% CI derived from PSA dataset (n = 10000)
intercepts <- paramDSA(4.167064, 3.526969, 4.784291)
log_scales <- paramDSA(0.321424, -0.06156224, 0.69504745)
GAMMAs <- 1/exp(log_scales)
lambdas <- (1/exp(intercepts))^(1/exp(log_scales))
n <- length(lambdas)

end <- 1200

# Calculate weibull survival
q <- matrix(sapply(c(1:n), function(x) exp(-lambdas[x]*((1:end)^GAMMAs[x]))), ncol = n)

# Calculate transition probabilities
p.1D <- matrix(NA, end, n)
p.1D[1, ] <- 1
for (i in 2:end){
  p.1D[i, ] <- ((q[i, ]))/(q[i-1, ])
}
1 - p.1D


# s <- matrix(NA, end,  15)
# s[1, ] <- 20000
# for (i in 2:end){
#   s[i, ] <- s[i-1, ] * (p.1D[i, ])
# }
# 
# par(mfrow = c(1, 1))
# plot(s[, 1]/s[1, 1], pch = 16, xlab = "Survival")
# par(new = TRUE)
# plot(s[, 2]/s[1, 2], pch = 16, axes = FALSE, ylab = NA, xlab = NA, col = "blue")
# par(new = TRUE)
# plot(s[, 15]/s[1, 3], pch = 16, axes = FALSE, ylab = NA, xlab = NA, col = "red")
# legend("right", horiz = FALSE,
#        legend = c("Base case", "Lower", "Upper"),
#        pch = c(16, 16, 16), col = c("black", "blue", "red")
# )
# apply(p.1Dis, 1, summary)
# apply(p.1Dis, 1, function (x) quantile(x, probs = c(0.025, 0.975)))

# par(mfrow = c(2, 1))
# apply(p.1Dis, 1, hist)
