p.1D    <- source(file.path("scripts", "p_1D_DSA.R"))[[1]]
p.10i <- as.matrix(read.csv2(file.path("PSA", "Output", "PSA", "P_10int_scale.csv")))[, -1]

# paramDSA <- function(x, min, max, steps = 7){
#   xout <- seq(min, max, length.out = steps*2 + 1)
#   xout <- xout[-(steps + 1)]
#   xout <- c(x, xout)
#   return(xout)
# }


lambdas <- paramDSA(0.8526, 0.8357, 0.8720)
n <- length(lambdas)
end <- 101

q <- sapply(c(1:n), function(x) exp(-lambdas[x]*1:2))



p.10temp <- q[2]/(q[1, ])

p.10 <- -4*log(1 - p.10temp)
p.10 <- 1 - exp(-p.10)

p.10 <- p.10 - (1 - p.1D)

p.10 

# s <- matrix(NA, 1200,  3)
# s[1, ] <- 20000
# for (i in 2:1200){
#   s[i, ] <- s[i-1, ] * (p.10[i, ])
# }
# 
# par(mfrow = c(1, 1))
# plot(s[, 1]/s[1, 1], pch = 16, xlab = "Survival")
# par(new = TRUE)
# plot(s[, 2]/s[1, 2], pch = 16, axes = FALSE, ylab = NA, xlab = NA, col = "blue")
# par(new = TRUE)
# plot(s[, 3]/s[1, 3], pch = 16, axes = FALSE, ylab = NA, xlab = NA, col = "red")
# legend("right", horiz = FALSE,
#        legend = c("Base case", "Lower", "Upper"),
#        pch = c(16, 16, 16), col = c("black", "blue", "red")
# )
summary(p.10i)
# 
# hist(p.10i)
