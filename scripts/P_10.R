# Inputdata for exponential
interceptpfs <- 2.922773
log_scalepfs <- 0.1892143
lambdapfs <- (1/exp(interceptpfs))^(1/exp(log_scalepfs))

# Exponential pfs
o <- exp(-lambdapfs*1:101)

# Calculate transition probabilities from survival
p.1Dpfs <- numeric(length(o))
p.1Dpfs[1] <- 0
for (i in 2:length(o)){
  p.1Dpfs[i] <- (o[i])/(o[i-1])
}

# Remove artificial first probability
p.1Dpfs <- p.1Dpfs[2:101]
p.1Dpfs <- 1 - p.1Dpfs


# Inputdata for weibull
interceptos <- 4.248068
log_scaleos <- 0.6365845
GAMMAos <- 1/exp(log_scaleos)
lambdaos <- (1/exp(interceptos))^(1/exp(log_scaleos))

# Weibull OS
oos <- exp(-lambdaos*((1:101)^GAMMAos))

# Calculate transition probabilities from survival
p.1Dos <- numeric(length(oos))
p.1Dos[1] <- 0
for (i in 2:length(oos)){
  p.1Dos[i] <- (oos[i])/(oos[i-1])
}

# Remove artificial first probability
p.1Dos <- p.1Dos[2:101]
p.1Dos <- 1 - p.1Dos

# Death probability == PFS - OS
p.10 <- p.1Dpfs - p.1Dos


# Convert yearly probabilities to monthly
monthrate2 <- -4*log(1-p.10)
monthtp2 <- 1-exp(-monthrate2)
p.10 <- rep(monthtp2, each = 12)
