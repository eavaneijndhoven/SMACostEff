# Read data
p.3D.DSAfun <- function(n.dsa){
  
  # paramDSA takes base, min, max value and amount of steps in sDSA, and returns extreme and intermediate values. 1: base case
  paramDSA <- function(x, min, max, steps = 7){
    xout <- seq(min, max, length.out = steps*2 + 1)
    xout <- xout[-(steps + 1)]
    xout <- c(x, xout)
    return(xout)
  }
  
  cbs <- read.csv2(file.path("Files", "Overlevingskansen__geslacht__leeftijd_28102019_115617.csv"))
  
  # 95% CI of rnorm(10000, 0, 0.001) with set.seed(29). This is also used in PSA
  mult <- paramDSA(0, -0.002001925, 0.001921450)
  
  # Convert percentages into proportions
  cbs$Survival <- cbs$Survival/100
  
  # Get mortality from survival
  cbs$Mortality <- 1 - cbs$Survival
  
  # Create PSA matrix
  # Since this gives a possibility of survival probability < 0 in early years, change those values to survival of previous year.
  cbs.DSA <- matrix(NA, 99, n.dsa)
  for(i in 1:n.dsa){
    cbs.DSA[, i] <- cbs$Mortality + mult[i]
  }
  for (i in 2:nrow(cbs.DSA)){
    cbs.DSA[i, ] <- ifelse(cbs.DSA[i, ] < 0, cbs.DSA[i - 1, ], cbs.DSA[i, ])
  }
  
  # Convert yearly survival probabilities to monthly probabilities
  cbs.mat <- matrix(NA, 1200, n.dsa)
  cbs.mat <- apply(cbs.DSA, 2, function(x) {x <- -1/12*log(1-x)
  x <- 1-exp(-x)
  x <- rep(x, each = 12)
  x <- append(x, x[1177:1188])}
  )
  
  return(cbs.mat)
}
################################################
## RUN FOR VISUAL INSPECTION OF DISTRIBUTION
################################################
# n.dsa <- 15
# 
# plot.DSA <- matrix(NA, 1200, n.dsa)
# plot.DSA[1, ] <- 21900
# 
# p.3Ddsa <- p.3D.DSAfun(n.dsa)
# 
# for (i in 2:nrow(plot.DSA)){
#   plot.DSA[i, ] <- plot.DSA[i - 1, ] * (1 - p.3Ddsa[i, ])
# }
# 
# palette <- rainbow(n.dsa)
# plot(plot.DSA[, 1]/21900)
# 
# for (i in 2:ncol(plot.DSA)){
#   par(new = TRUE)
#   plot(plot.DSA[, i]/21900, axes = FALSE, ylab = NA, xlab = NA, col = palette[i])
# }


