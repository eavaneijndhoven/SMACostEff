# Inputdata for weibull
intercept <- 4.167065
log_scale <- 0.321424
GAMMA <- 1/exp(log_scale)
lambda <- (1/exp(intercept))^(1/exp(log_scale))

# Weibull survival
x <- exp(-lambda*((1:(n.t + 1))^GAMMA))


# Calculate transition probabilities from survival
p.1D <- numeric(length(x))
p.1D[1] <- 0
for (i in 2:length(x)){
  p.1D[i] <- (x[i])/(x[i-1])
}
# Remove artificial first probability
p.1D <- p.1D[2:(n.t + 1)]
p.1D <- 1 - p.1D