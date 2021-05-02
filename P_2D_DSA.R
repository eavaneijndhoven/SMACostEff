# p.2Dis <- as.matrix(read.csv2(file.path("PSA", "Output", "PSA", "P_2Dint_scale.csv")))[, -1]
# dimnames(p.2Dis) <- list(c("Intercept", "Log_scale"))
# paramDSA <- function(x, min, max, steps = 7){
#   xout <- seq(min, max, length.out = steps*2 + 1)
#   xout <- xout[-(steps + 1)]
#   xout <- c(x, xout)
#   return(xout)
# }

# Range of intercepts and log_scales: 95% CI derived from PSA dataset (n = 10000)
intercepts <- paramDSA(3.833436, 3.645634, 4.014528) 
log_scales <- paramDSA(-0.6874054, -0.9068469, -0.4779306)
GAMMAs <- 1/exp(log_scales)
lambdas <- (1/exp(intercepts))^(1/exp(log_scales))
n <- length(lambdas)
end <- 100

# Calculate weibull survival
q <- matrix(sapply(c(1:n), function(x) exp(-lambdas[x]*((1:(end + 1))^GAMMAs[x]))), ncol = n)

# Calculate transition probabilities from survival
p.2D <- matrix(NA, end + 1, n)
p.2D[1, ] <- 1
for (i in 2:(end + 1)){
p.2D[i, ] <- q[i, ]/q[i-1, ]
}

# Remove artificial first probability
p.2D <- p.2D[2:(end + 1), ]
p.2D <- 1-p.2D

p.2D.mat <- matrix(NA, 1200, ncol = n)
# Convert yearly probabilities to monthly
for (i in 1:n){  
  p.2D[, i] <- -1/12*log(1-p.2D[, i])
  p.2D[, i] <- 1-exp(-p.2D[, i])
  p.2D.mat[, i] <- rep(p.2D[, i], each = 12)
}

p.2D.mat

# 
# s <- matrix(NA, 1200,  15)
# s[1, ] <- 20000
# for (i in 2:1200){
#   s[i, ] <- s[i-1, ] * (1 - p.2D.mat[i, ])
# }
# 
# par(mfrow = c(1, 1))
# plot(s[, 1]/s[1, 1], pch = 16, xlab = "Survival")
# par(new = TRUE)
# plot(s[, 2]/s[1, 2], pch = 16, axes = FALSE, ylab = NA, xlab = NA, col = "blue")
# par(new = TRUE)
# plot(s[, 3]/s[1, 15], pch = 16, axes = FALSE, ylab = NA, xlab = NA, col = "red")
# legend("right", horiz = FALSE,
#        legend = c("Base case", "Lower", "Upper"),
#        pch = c(16, 16, 16), col = c("black", "blue", "red")
# )
# apply(p.2Dis, 1, summary)
# apply(p.2Dis, 1, function (x) quantile(x, probs = c(0.025, 0.975)))
# 
# par(mfrow = c(2, 1))
# apply(p.2Dis, 1, hist)
