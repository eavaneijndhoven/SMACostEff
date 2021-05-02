# Inputdata for weibull
intercept <- 4.248068
log_scale <- 0.6365845
GAMMA <- 1/exp(log_scale)
lambda <- (1/exp(intercept))^(1/exp(log_scale))

# Weibull survival
t <- exp(-lambda*((1:101)^GAMMA))

# Calculate transition probabilities
p.1D <- numeric(length(t))
p.1D[1] <- 0
for (i in 2:length(t)){
  p.1D[i] <- (t[i])/(t[i-1])
}

# Remove artificial first probability
p.1D <- p.1D[2:101]
p.1D <- 1 - p.1D

# Convert yearly probabilities to monthly
monthrate2 <- -4*log(1-p.1D)
monthtp2 <- 1-exp(-monthrate2)
p.1D <- rep(monthtp2, each = 12)
