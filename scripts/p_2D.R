# Inputdata for weibull
intercept <- 3.833436
log_scale <- -0.6874054
GAMMA <- 1/exp(log_scale)
lambda <- (1/exp(intercept))^(1/exp(log_scale))

# Weibull survival
q <- exp(-lambda*((1:((n.t/12)+1))^GAMMA))

# Calculate transition probabilities from survival
p.2D <- numeric(length(q))
p.2D[1] <- 1
for (i in 2:length(q)){
  p.2D[i] <- q[i]/q[i-1]
}

# Remove artificial first probability
p.2D <- p.2D[2:((n.t/12)+1)]
p.2D <- 1-p.2D

# Convert yearly probabilities to monthly
monthrate2 <- -1/12*log(1-p.2D)
monthtp2 <- 1-exp(-monthrate2)
p.2D <- rep(monthtp2, each = 12)