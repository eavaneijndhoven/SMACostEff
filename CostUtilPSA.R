# Random draws for health state costs
CostPSA <- function(n.psa, Trt){
  # Base case costs depending on treatment setting
  ifelse(Trt == TRUE, psa.state <- list(psa.1 = 9625*1.0125*1.0196, psa.2 = 10197*1.0125*1.0196,
                      psa.3 = 5680.5*1.0125*1.0196, psa.0 = 14725.7*1.0125*1.0196),
                      psa.state <- list(psa.1 = 10701*1.0125*1.0196, psa.2 = 10457*1.0125*1.0196,
                      psa.3 = 5808*1.0125*1.0196, psa.0 = 14725.7*1.0125*1.0196))
  # Calculate SE
  psa.state <- lapply(psa.state, function (x) {
    append(x, (x*1.25-x*0.75)/(2*1.96))
    
  })
  
  # Calculate gamma parameters
  psa.state <- lapply(psa.state, function(x){
    append(x, c((x[1]/x[2])^2, x[2]^2/x[1]))
  })
  # Randomize with gamma distribution
  psa.c <- vector("list", 4)
  psa.c <- lapply(psa.state, function(x){
    x <- rgamma(n.psa, x[3], 1/x[4])
  })
  
  return(psa.c)
}




UtilPSA <- function(n.psa){
  # Input basecase values, with SEM derived from range in ZIN report
  psa.state <- list(psa.1 = c(0.733, 0.01 ), psa.2 = c(0.752, 0.015), psa.3 = c(0.878, 0.028), psa.0 = c(0.733, 0.001))
  
  # Calculate alpha
  psa.state <- lapply(psa.state, function(x){
    append(x, x[1]^2 * (1 - x[1])/(x[2]^2) - 1 * x[1])
    
  })
  
  # Calculate beta
  psa.state <- lapply(psa.state, function(x){
    append(x, x[3] * (1 - x[1]) / x[1])
  })
  
  # Create random draws based on parameters
  psa.u <- vector("list", 4)
  psa.u <- lapply(psa.state, function(x){
    x <- rbeta(n.psa, x[3], x[4])
  })
  
  return(psa.u)
}

# Randomize treatment costs with gamma distribution
TrtCost <- function(n.psa, Treat = c("Spinraza", "Zolgensma")){
  ifelse(Treat == "Spinraza", base <- 83300*1.0125*1.0196, base <- 2000000)
  Trtmin   <- base * 0.75
  Trtmax   <- base * 1.25
  Trtsem   <- (Trtmax - Trtmin) / (1.96^2)
  Trtalpha <- (base/Trtsem)^2
  Trtgamma <- (Trtsem^2)/base
  
  Costout <- rgamma(n.psa, Trtalpha, 1/Trtgamma)
  return(Costout)
}

### RUN TO CHECK DISTRIBUTION

# foo <- UtilPSA(10000)
# 
# state <- c("1", '2', "3", "0")
# 
# par(mfrow = c(2,2))
# for (i in 1:4){
#   hist(foo[[i]], main = paste0("Histogram of u.", state[i]), xlab = paste0("u.", state[i]))
# }
# 
# 
# bar <- CostPSA(10000)
# 
# state <- c("1", '2', "3", "0")
# 
# par(mfrow = c(2,2))
# for (i in 1:4){
#   hist(bar[[i]], main = paste0("Histogram of c.", state[i]), xlab = paste0("c.", state[i]))
# }
# 
# tcost <- ZolgCost(10000)
# hist(tcost)
