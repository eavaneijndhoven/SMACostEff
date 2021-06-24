# RANDOMIZATION OF TRANSITION PROBABILITIES FOR PROBABILISTIC SENSITIVITY ANALYSIS

# Random draws for p.1D PSA
p.1D.PSAfun <- function(n.psa){
  
  # Cholesky matrix with intercept and log_scale
  cholmat <- matrix(c(0.32062961, 0.00000000, 0.05216662, 0.1832504),
                    2, 2, byrow = TRUE)
  intercept <- 4.16706450522185
  log_scale <- 0.321424017677827
  
  # Create empty matrix with row 1 = intercept and row 2 = log_scale. Each column represents input for 1 Monte Carlo simulation
  rands <- matrix(NA, nrow = 2, ncol = n.psa)
  
  # Fill matrix with randomly generated intercept and log_scale
  rands <- apply(rands, 2, function (x) c(intercept, log_scale) + cholmat %*% rnorm(2, mean = 0, sd = 1))
  
  # Create empty matrix with row 1 = gamma and row 2 = lambda. Each column represents input for 1 Monte Carlo simulation
  rands.wei <- matrix(NA, nrow = 2, ncol = n.psa)
  rands.wei[1, ] <- 1/exp(rands[2, ])
  rands.wei[2, ] <- (1/exp(rands[1, ]))^(1/exp(rands[2, ]))
  
  # Create matrix for survival probabilities
  q.PSA <- matrix(NA, nrow = 1201, ncol = n.psa)
  
  #Calculate survival using rands.wei
  for (i in 1:n.psa){
   q.PSA[, i] <- exp(-rands.wei[2, i]*((1:1201)^rands.wei[1, i]))
  }

  # Calculate transition probabilities from survival
  p.1D.PSA <- matrix(NA, nrow = 1201, ncol = n.psa)
  
  for (i in 1:n.psa){  
    p.1D.PSA[1, i] <- 1
    for (j in 2:length(p.1D.PSA[, i])){
      p.1D.PSA[j, i] <- ((q.PSA[j, i]))/(q.PSA[j-1, i])
    }
  }
  
  # Remove artificial first probability, and convert to death probability
  p.1D.PSA <- p.1D.PSA[2:1201, ]
  p.1D.PSA <- 1-p.1D.PSA
  
  return(p.1D.PSA)
}




# Random draws for p.2D PSA
p.2D.PSAfun <- function(n.psa){
  
  # Cholesky matrix with intercept and log_scale
  cholmat <- matrix(c(0.09407148, 0.00000000, 0.05556044, 0.09199783),
                    2, 2, byrow = TRUE)
  intercept <- 3.833436
  log_scale <- -0.6874054
  
  # Create empty matrix with row 1 = intercept and row 2 = log_scale. Each column represents input for 1 Monte Carlo simulation
  rands <- matrix(NA, nrow = 2, ncol = n.psa)
  
  # Fill matrix with randomly generated intercept and log_scale
  rands <- apply(rands, 2, function (x) c(intercept, log_scale) + cholmat %*% rnorm(2, mean = 0, sd = 1))
  
  # Create empty matrix with row 1 = gamma and row 2 = lambda. Each column represents input for 1 Monte Carlo simulation
  rands.wei <- matrix(NA, nrow = 2, ncol = n.psa)
  rands.wei[1, ] <- 1/exp(rands[2, ])
  rands.wei[2, ] <- (1/exp(rands[1, ]))^(1/exp(rands[2, ]))
  
  # Create matrix for survival probabilities
  q.PSA <- matrix(NA, nrow = 101, ncol = n.psa)
  
  # Calculate survival using rands.wei
  for (i in 1:n.psa){
    q.PSA[, i] <- exp(-rands.wei[2, i]*((1:101)^rands.wei[1, i]))
  }
  # Calculate transition probabilities from survival
  p.2D.PSA <- matrix(NA, nrow = 101, ncol = n.psa)
  
  for (i in 1:n.psa){  
    p.2D.PSA[1, i] <- 1
    for (j in 2:length(p.2D.PSA[, i])){
      p.2D.PSA[j, i] <- ((q.PSA[j, i]))/(q.PSA[j-1, i])
    }
  }
  
  # Remove artificial first probability, and convert to death probability
  p.2D.PSA <- p.2D.PSA[2:101, ]
  p.2D.PSA <- 1-p.2D.PSA
  
  p.2D.mat <- matrix(NA, 1200, ncol = n.psa)
  # Convert yearly probabilities to monthly
  for (i in 1:n.psa){  
    p.2D.PSA[, i] <- -1/12*log(1-p.2D.PSA[, i])
    p.2D.PSA[, i] <- 1-exp(-p.2D.PSA[, i])
    p.2D.mat[, i] <- rep(p.2D.PSA[, i], each = 12)
  }
  
  return(p.2D.mat)
}


# Random draws for p.3D PSA
p.3D.PSAfun <- function(n.psa){
  
  cbs <- read.csv2(file.path("Files", "Overlevingskansen__geslacht__leeftijd_28102019_115617.csv"))
  
  
  # Convert percentages into proportions
  cbs$Survival <- cbs$Survival/100
  
  # Get mortality from survival
  cbs$Mortality <- 1 - cbs$Survival
  
  # Create PSA matrix, and add normally distributed noise with mean = 0 and sd = 0.01 to survival.
  # Since this gives a possibility of survival probability < 0 in early years, change those values to survival of previous year.
  cbs.PSA <- matrix(NA, 99, n.psa)
  cbs.PSA <- apply(cbs.PSA, 2, function(x) cbs$Mortality + rnorm(1, 0, 0.001))
  for (i in 2:nrow(cbs.PSA)){
    cbs.PSA[i, ] <- ifelse(cbs.PSA[i, ] < 0, cbs.PSA[i - 1, ], cbs.PSA[i, ])
  }
  
  # Convert yearly survival probabilities to monthly probabilities
  cbs.mat <- matrix(NA, 1200, n.psa)
  cbs.mat <- apply(cbs.PSA, 2, function(x) {x <- -1/12*log(1-x)
  x <- 1-exp(-x)
  x <- rep(x, each = 12)
  x <- append(x, x[1177:1188])}
  )
  
  return(cbs.mat)
}

# Random draws for p.1D PSA in BSC arm
p.1D.BSC.PSAfun <- function(n.psa){
  
  # Cholesky matrix with intercept and log_scale
  cholmat <- matrix(c(0.5446921, 0.00000000, 0.1099314, 0.1916723),
                    2, 2, byrow = TRUE)
  intercept <- 4.248068
  log_scale <- -0.6365845
  
  # Create empty matrix with row 1 = intercept and row 2 = log_scale. Each column represents input for 1 Monte Carlo simulation
  rands <- matrix(NA, nrow = 2, ncol = n.psa)
  
  # Fill matrix with randomly generated intercept and log_scale
  rands <- apply(rands, 2, function (x) c(intercept, log_scale) + cholmat %*% rnorm(2, mean = 0, sd = 1))
  
  # Create empty matrix with row 1 = gamma and row 2 = lambda. Each column represents input for 1 Monte Carlo simulation
  rands.wei <- matrix(NA, nrow = 2, ncol = n.psa)
  rands.wei[1, ] <- 1/exp(rands[2, ])
  rands.wei[2, ] <- (1/exp(rands[1, ]))^(1/exp(rands[2, ]))
  
  # Create matrix for survival probabilities
  q.PSA <- matrix(NA, nrow = 101, ncol = n.psa)
  
  # Calculate survival using rands.wei
  for (i in 1:n.psa){
    q.PSA[, i] <- exp(-rands.wei[2, i]*((1:101)^rands.wei[1, i]))
  }
  # Calculate transition probabilities from survival
  p.1D.BSC.PSA <- matrix(NA, nrow = 101, ncol = n.psa)
  
  for (i in 1:n.psa){  
    p.1D.BSC.PSA[1, i] <- 1
    for (j in 2:length(p.1D.BSC.PSA[, i])){
      p.1D.BSC.PSA[j, i] <- ((q.PSA[j, i]))/(q.PSA[j-1, i])
    }
  }
  
  # Remove artificial first probability, and convert to death probability
  p.1D.BSC.PSA <- p.1D.BSC.PSA[2:101, ]
  p.1D.BSC.PSA <- 1-p.1D.BSC.PSA
  
  p.1D.BSC.mat <- matrix(NA, 1200, ncol = n.psa)
  # Convert yearly probabilities to monthly
  for (i in 1:n.psa){  
    p.1D.BSC.PSA[, i] <- -4*log(1-p.1D.BSC.PSA[, i])
    p.1D.BSC.PSA[, i] <- 1-exp(-p.1D.BSC.PSA[, i])
    p.1D.BSC.mat[, i] <- rep(p.1D.BSC.PSA[, i], each = 12)
  }
  
  return(p.1D.BSC.mat)
}

# Random draws for p.10 PSA in BSC arm
p.10.BSC.PSAfun <- function(n.psa, p.1Dbsc){
  
 
  # Cholesky matrix with intercept
  cholmat <- matrix(c(0.1892143),
                    1, 1, byrow = TRUE)
  intercept <- 2.922773
  
  # Create empty matrix with row 1 = intercept. Each column represents input for 1 Monte Carlo simulation
  rands <- matrix(NA, nrow = 1, ncol = n.psa)
  
  # Fill matrix with randomly generated intercept
  rands <- apply(rands, 2, function (x) intercept + cholmat %*% rnorm(2, mean = 0, sd = 1))
  
  # Create empty matrix with row 1 = gamma and row 2 = lambda. Each column represents input for 1 Monte Carlo simulation
  rands.wei <- matrix(NA, nrow = 1, ncol = n.psa)
  rands.wei[1, ] <- (1/exp(rands[1, ]))^(1/exp(rands[2, ]))
  
  # Create matrix for survival probabilities
  q.PSA <- matrix(NA, nrow = 101, ncol = n.psa)
  
  # Calculate survival using rands.wei
  for (i in 1:n.psa){
    q.PSA[, i] <- exp(-rands.wei[1, i]*1:101)
  }
  # Calculate transition probabilities from survival
  p.1DPFS.BSC.PSA <- matrix(NA, nrow = 101, ncol = n.psa)
  
  for (i in 1:n.psa){  
    p.1DPFS.BSC.PSA[1, i] <- 1
    for (j in 2:length(p.1DPFS.BSC.PSA[, i])){
      p.1DPFS.BSC.PSA[j, i] <- ((q.PSA[j, i]))/(q.PSA[j-1, i])
    }
  }
  
  # Remove artificial first probability, and convert to death probability
  p.1DPFS.BSC.PSA <- p.1DPFS.BSC.PSA[2:101, ]
  p.1DPFS.BSC.PSA <- 1-p.1DPFS.BSC.PSA
  
  
  p.1DPFS.BSC.mat <- matrix(NA, 1200, ncol = n.psa)
  # Convert yearly probabilities to monthly
  for (i in 1:n.psa){  
    p.1DPFS.BSC.PSA[, i] <- -4*log(1-p.1DPFS.BSC.PSA[, i])
    p.1DPFS.BSC.PSA[, i] <- 1-exp(-p.1DPFS.BSC.PSA[, i])
    p.1DPFS.BSC.mat[, i] <- rep(p.1DPFS.BSC.PSA[, i], each = 12)
  }
  
  p.10.PSA <- p.1DPFS.BSC.mat - p.1Dbsc
  
  
  
  return(p.10.PSA)
}
