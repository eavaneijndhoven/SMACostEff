packs <- c("foreach", "doSNOW", "devtools", "ggplot2")
for(i in 1:length(packs)){
  if (!packs[i] %in% installed.packages()) {install.packages(packs[i])}
}
library(foreach)
library(doSNOW)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyverse)
source(file.path("scripts", "gg_pubthemes.R"))



##################################### Model input ##########################################

# Model input
n.i   <- 10000 # number of simulated individuals
n.t   <- 1200                               # time horizon, 1200 cycles
v.n   <- c("1","2","3", "0",                   # Health states in the model
           "Death")
n.s   <- length(v.n)                           # the number of health states
set.seed(29)

# Generate DSA inputs based on min, max and number of steps
paramDSA <- function(x, min, max, steps = 7){
  xout <- seq(min, max, length.out = steps*2 + 1)
  xout <- xout[-(steps + 1)]
  xout <- c(x, xout)
  return(xout)
}


TypeRange <- "95% C.I."
# Input survival
p.2D    <- source(file.path("scripts", "p_2D_DSA.R"))[[1]]
source(file.path("scripts", "P_3D_DSA.R"))
p.3D    <- p.3D.DSAfun(15)
p.0D    <- source(file.path("scripts", "p_1D_DSA.R"))[[1]]
p.10bsc <- source(file.path("scripts", "p_10_DSA.R"))[[1]]
p.10    <- matrix(rep(0, n.t*n), ncol = n)
p.1D    <- source(file.path("scripts", "p_1D_DSA.R"))[[1]]
p.1Dbsc <- source(file.path("scripts", "p_1D_BSC_DSA.R"))[[1]]
p.2D    <- ifelse(p.2D > p.3D, p.2D, p.3D)


# Input costs and utilities
c.1     <- paramDSA(9625, 6045, 10739)*1.0125*1.0196
c.2     <- paramDSA(10197, 7647, 12746)*1.0125*1.0196
c.3     <- paramDSA(5680, 4260, 7100)*1.0125*1.0196
c.0     <- paramDSA(14725, 13882, 16490)*1.0125*1.0196
c.1bsc  <- paramDSA(10701, 10053, 11948)*1.0125*1.0196

c.zolg  <- paramDSA(2000000, 1500000, 2500000)
c.spin  <- paramDSA(83300, 62475, 104125)

u.1     <- paramDSA(0.733, 0.714, 0.753)
u.2     <- paramDSA(0.752, 0.721, 0.783)
u.3     <- paramDSA(0.878, 0.821, 0.993)
u.0     <- paramDSA(0.733, 0.714, 0.753)

# Discount
d.c   <- paramDSA(0.04, 0, 0.08)                                  # discounting of costs by 4%
d.e   <- paramDSA(0.015, 0, 0.03)                                 # discounting of QALYs by 1.5%

# Create discount matrix
m.dwc <- sapply(d.c, function (x) c(1, rep(1 / (1 + x) ^ (0:((n.t/12) -1)), each = 12)))    # calculate the cost discount weight based on the discount rate d.c 
m.dwe <- sapply(d.e, function (x) c(1, rep(1 / (1 + x) ^ (0:((n.t/12) -1)), each = 12)))    # calculate the QALY discount weight based on the discount rate d.e

#################################### Create DSA input parameters ##########################################

# dsa.rownames: names of parameters included in PSA
# dsa.labels: deviations from base case values
# dsa.levels: number of labels
# dsa.pos: variable to place dsa.labels in DSA matrix
# n.dsa: umber of input variables for DSA
# DSA: matrix containing information which input parameters to use in each run of the DSA. In each run, one of the 
# adjusted variables is used, the rest are kept at base case level. Each column contains the information for one DSA run.
# n.dsasims: total number of runs needed to complete the DSA. 


dsa.rownames <- c("d.c", "d.e", "c.zolg", "c.spin", "c.1", "c.2", "c.3", "c.0", "u.1", "u.2", "u.3", "u.0", "p.1D", "p.2D", "p.3D", "p.0D", "p.10")
dsa.fullnames <- c("Discounting of\ncosts", "Discounting of\neffects", "Price Zolgensma", "Price Spinraza", "Cost state 1",
                   "Cost state 2", "Cost state 3", "Cost perm.\nventilation", "Utility state 1",
                   "Utility state 2", "Utility state 3", "Utility perm.\nventilation",
                   "Survival state 1", "Survival state 2", "Survival state 3", "Survival perm.\nventilation",
                   "Prob. state 1 to\nperm. ventilation")
dsa.labels <- c("lower", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "upper")
dsa.levels <- length(dsa.labels)
dsa.pos    <- 2:(dsa.levels + 1)
dsa.colnames <- c("base", paste(rep(dsa.rownames, each = dsa.levels), dsa.labels, sep = "_"))
n.dsa <- length(dsa.rownames) 
DSA <- matrix(1, nrow = n.dsa, ncol = (n.dsa * dsa.levels) + 1, 
              dimnames = list(dsa.rownames, dsa.colnames))
n.dsasims <- length(dsa.colnames)
for(r in dsa.pos){
  for(i in seq(r, n.dsasims, dsa.levels)){
    ifelse(r < dsa.levels, DSA[(i/dsa.levels) + 1, i] <- r, DSA[i/dsa.levels, i] <- r)
  }
}



##################################### Functions ###########################################

# THE NEW samplev() FUNCTION
# efficient implementation of the rMultinom() function of the Hmisc package #### 

samplev <- function (probs, m) {  
  d <- dim(probs)
  n <- d[1]
  k <- d[2]
  lev <- dimnames(probs)[[2]]
  if (!length(lev)) 
    lev <- 1:k
  ran <- matrix(lev[1], ncol = m, nrow = n)
  U <- t(probs)
  for(i in 2:k) {
    U[i, ] <- U[i, ] + U[i - 1, ]
  }
  if (any((U[k, ] - 1) > 1e-05))
    stop("error in multinom: probabilities do not sum to 1")
  
  for (j in 1:m) {
    un <- rep(runif(n), rep(k, n))
    ran[, j] <- lev[1 + colSums(un > U)]
  }
  ran
}
# The MicroSim function keeps track of what happens to each individual during each cycle. 

MicroSimDSA <- function(n.i, n.t, v.n, f, c.Trt,
                        treat = c("Zolgensma", "BSC", "Spinraza")) {
  # Arguments:  
  # v.M_1:    vector of initial states for individuals
  # n.i:      number of individuals
  # n.t:      total number of cycles to run the model
  # v.n:      vector of health state names 
  # d.c:      discount rate for costs
  # d.e:      discount rate for health outcome (QALYs)
  # TR.out:   should the output include a Microsimulation trace? (default is TRUE). Is used for trace plots with simplot()
  # TS.out:   should the output include a matrix of transitions between states? (default is TRUE)
  # data.out: should the output include the matrices with states, costs and effects (default is TRUE) (adds +1 GB data in this model, only run when needed)
  # treat:    which treatment setting should be used for the model? Choice between Zolgensma and bscraza
  # seed:     starting seed number for random number generator (default is 1) 
  # Makes use of:
  # Probs:   function for the estimation of transition probabilities
  # Costs:   function for the estimation of cost state values
  # Effs:    function for the estimation of state specific health outcomes (QALYs)
  
  if (treat == "Zolgensma") {v.M_1 <- c(rep("1", n.i*0.08), rep("2", n.i*0.59), rep("3", n.i*0.33)); Trt <- TRUE}
  if (treat == "BSC")       {v.M_1 <- c(rep("1", n.i*0.68), rep("0", n.i*0.32)); p.10 <- p.10bsc; c.1 <- c.1bsc; p.1D <- p.1Dbsc; Trt <- FALSE}
  if (treat == "Spinraza")  {v.M_1 <- c(rep("1", n.i*0.5641), rep("2", n.i*0.1859), rep("3", n.i*0), rep("0", n.i*0.25)); Trt <- TRUE}
  
  
  
  # Create the matrix capturing the state name/costs/health outcomes for all individuals at each time point 
  m.M <- m.C <- m.E <-  matrix(nrow = n.i, ncol = n.t + 1, 
                               dimnames = list(paste("ind", 1:n.i, sep = " "), 
                                               paste("cycle", 0:n.t, sep = " ")))  
  
  m.M[, 1] <- v.M_1                     # indicate the initial health state   
  
  m.C[, 1] <- Costs(m.M[, 1], Trt, t = 1, prev = NA, c.Trt)  # estimate costs per individual for the initial health state
  if (treat == "Zolgensma") {
    m.C[, 1] <- m.C[, 1] + c.zolg[DSA[3, f]] + 1336.11
  }  # add treatment costs of Zolgensma   
  if (treat == "Spinraza"){
    m.C[, 1] <- m.C[, 1] + (c.Trt[DSA[4, f]]*4) 
  }
  
  m.E[, 1] <- Effs(m.M[, 1], f = f)  # estimate QALYs per individual for the initial health state
  
  
  
  for (t in 1:n.t) {
    m.p <- Probs(m.M[, t], t, f)           # calculate the transition probabilities at cycle t 
    m.M[, t + 1] <- samplev(prob = m.p, m = 1)  # sample the next health state and store that state in matrix m.M 
    m.C[, t + 1] <- Costs(m.M[, t + 1], Trt, t = t, prev = m.M[, t], c.Trt, f = f)   # estimate costs per individual during cycle t + 1 conditional on treatment
    m.E[, t + 1] <- Effs(m.M[, t + 1], f = f)   # estimate QALYs per individual during cycle t + 1 conditional on treatment
    cat('\r', paste(round(t/n.t * 100, 2), "% done    ", sep = " "))       # display the progress of the simulation
    
  } # close the loop for the time points 
  
  
  tc <- m.C %*% m.dwc[, DSA[1, f]]            # total (discounted) cost per individual
  te <- (m.E/12) %*% m.dwe[, DSA[2, f]]       # total (discounted) QALYs per individual 
  
  tc_hat <- mean(tc)       # average (discounted) cost 
  te_hat <- mean(te)        # average (discounted) QALYs
  
  results <- list(tc_hat = tc_hat, te_hat = te_hat) # store the results from the simulation in a list  
  return(results)  # return the results
  
  
}  # end of the MicroSim function  


#### Probability function
# The Probs function that updates the transition probabilities of every cycle is shown below.

Probs <- function(M_it, t = t, f = f) { 
  # M_it:    health state occupied by individual i at cycle t (character variable)
  
  m.p.it <- matrix(NA, n.s, n.i)     # create vector of state transition probabilities
  rownames(m.p.it) <- v.n  # assign vector names
  
  # update the v.p with the appropriate probabilities
  m.p.it[,M_it == "1"]     <- c(1 - (p.1D[t,  DSA[13, f]] + p.10[t, DSA[17, f]]), 0, 0, p.10[t, DSA[17, f]], p.1D[t, DSA[13, f]])                 
  m.p.it[,M_it == "2"]     <- c(0, (1 - p.2D[t, DSA[14, f]]), 0, 0, p.2D[t, DSA[14, f]])   
  m.p.it[,M_it == "3"]     <- c(0, 0, 1 - p.3D[t, DSA[15, f]], 0, p.3D[t, DSA[15, f]])                           
  m.p.it[,M_it == "0"]     <- c(0, 0, 0, 1 - p.0D[t, DSA[16, f]], p.0D[t, DSA[16, f]])
  m.p.it[,M_it == "Death"] <- c(0, 0, 0, 0, 1)
  
  
  ifelse(colSums(m.p.it) == 1, return(t(m.p.it)), print("Probabilities do not sum to 1"))    # return the transition probabilities or produce an error
  
}


### Costs function
# The Costs function estimates the costs at every cycle.

Costs <- function (M_it, Trt, t, prev, c.Trt, f) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  
  
  c.it <- 0                                  # by default the cost for everyone is 0
  c.it[M_it == "1"]     <- ifelse(t %in% seq(5, 1200, by = 4),  c.1[DSA[5, f]] + c.Trt[DSA[4, f]] * Trt, c.1[DSA[5, f]])               
  c.it[M_it == "2"]     <- ifelse(t %in% seq(5, 1200, by = 4),  c.2[DSA[6, f]] + c.Trt[DSA[4, f]] * Trt, c.2[DSA[6, f]])
  c.it[M_it == "3"]     <- ifelse(t %in% seq(5, 1200, by = 4),  c.3[DSA[7, f]] + c.Trt[DSA[4, f]] * Trt, c.3[DSA[7, f]])
  c.it[M_it == "0"]     <- c.0[DSA[8, f]]                
  c.it[M_it == "Death"] <- 0
  c.it[M_it == "Death" & prev == "1"] <- c.0[DSA[8, f]]*3    # dying gives a one time cost
  
  return(c.it)        		                   # return the costs
  
}


### Health outcome function 
# The Effs function to update the utilities at every cycle.

Effs <- function (M_it, cl = 1, f) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  # Trt:  is the individual treated? (default is FALSE) 
  # cl:   cycle length (default is 1)
  
  u.it <- 0                      # by default the utility for everyone is zero
  u.it[M_it == "1"]     <- u.1[DSA[9, f]]      
  u.it[M_it == "2"]     <- u.2[DSA[10, f]] 
  u.it[M_it == "3"]     <- u.3[DSA[11, f]]
  u.it[M_it == "0"]     <- u.0[DSA[12, f]]      
  u.it[M_it == "Death"] <- 0
  
  QALYs <-  u.it * cl            # calculate the QALYs during cycle t
  return(QALYs)                  # return the QALYs
  
}


##################################### Run the simulation ##################################
# Setup multicore cluster
cores <- parallel::detectCores()
cl <- makeSOCKcluster(cores[1] - 1)
registerDoSNOW(cl)

# START SIMULATION
zolgout <- matrix(NA, n.dsasims, 2)
zolgtemp <- numeric(2)
bscout  <- matrix(NA, n.dsasims, 2)
bsctemp <- numeric(2)
spinout <- matrix(NA, n.dsasims, 2)
spintemp <- numeric(2)


n.dsasims <- 4
# Setup progress bar
pb <- txtProgressBar(min = 1, max = n.dsasims, style = 3)
progress <- function(n) {setTxtProgressBar(pb, n)}
opts <- list(progress = progress)

p <- Sys.time()
p

zolgout <- foreach(f = 1:n.dsasims, .combine = rbind, .options.snow = opts) %dopar% {
  sim_zolg  <- MicroSimDSA(n.i, n.t, v.n, f = f, treat = "Zolgensma",
                           c.Trt = rep(0, dsa.levels + 1)) # run for Zolgensma
  
  zolgtemp[1] <- sim_zolg$tc_hat
  zolgtemp[2] <- sim_zolg$te_hat
  
  setTxtProgressBar(pb, f)
  zolgtemp
}
close(pb)

pb <- txtProgressBar(min = 1, max = n.dsasims, style = 3)
spinout <- foreach(f = 1:n.dsasims, .combine = rbind, .options.snow = opts) %dopar% {
  sim_spin <- MicroSimDSA(n.i, n.t, v.n, f = f, treat = "Spinraza",
                          c.Trt = c.spin) # Run for Spinraza
  
  spintemp[1] <- sim_spin$tc_hat
  spintemp[2] <- sim_spin$te_hat
  
  setTxtProgressBar(pb, f)
  spintemp
}
close(pb)

pb <- txtProgressBar(min = 1, max = n.dsasims, style = 3)
bscout <- foreach(f = 1:n.dsasims, .combine = rbind, .options.snow = opts) %dopar% {
  sim_bsc  <- MicroSimDSA(n.i, n.t, v.n, f = f, treat = "BSC",
                          c.Trt = rep(0, dsa.levels + 1)) # run for BSC
  
  bsctemp[1]  <- sim_bsc$tc_hat
  bsctemp[2]  <- sim_bsc$te_hat
  
  setTxtProgressBar(pb, f)
  bsctemp
}
close(pb)


stopCluster(cl)
comp.time <- Sys.time() - p
comp.time

# Combine results of DSAs to calculate incremental costs and effects and ICERs
zolgbsccost <- c(zolgout[, 1] - bscout[, 1])
zolgbsceff  <- c(zolgout[, 2] - bscout[, 2])
zolgbscICER <- zolgbsccost / zolgbsceff

spinbsccost <- c(spinout[, 1] - bscout[, 1])
spinbsceff  <- c(spinout[, 2] - bscout[, 2])
spinbscICER <- spinbsccost / spinbsceff

zolgspincost <- c(zolgout[, 1] - spinout[, 1])
zolgspineff  <- c(zolgout[, 2] - spinout[, 2])
zolgspinICER <- zolgspincost / zolgspineff


# Create labels for sDSA dataframe
Labels <- vector("list", n.dsa)
labeltemp <- rep(NA, dsa.levels)
for (i in 1:n.dsa){
  labeltemp <- c(as.character(min(get(dsa.rownames[i]))), rep(NA, dsa.levels - 2), as.character(max(get(dsa.rownames[i]))))
  Labels[[i]] <- labeltemp
}
Labels <- unlist(Labels)

# Data.frames for sDSA script
zolgbsc <- data.frame(Scenario = 1:n.dsasims,
                      Parameter = c("base case", rep(dsa.rownames, each = dsa.levels)),
                      ScenarioNumb = c(1, rep(2:n, n.dsa)),
                      Bound = c("base case", rep(c(rep("lower", (dsa.levels)/2), rep("upper", dsa.levels/2)), n.dsa)),
                      ICER = zolgbscICER,
                      Labels = c("base case", Labels),
                      TypeRange = TypeRange,
                      FullName = c("Base case", rep(c(dsa.fullnames), each = dsa.levels)))

spinbsc <- data.frame(Scenario = 1:n.dsasims,
                      Parameter = c("base case", rep(dsa.rownames, each = dsa.levels)),
                      ScenarioNumb = c(1, rep(2:n, n.dsa)),
                      Bound = c("base case", rep(c(rep("lower", (dsa.levels)/2), rep("upper", dsa.levels/2)), n.dsa)),
                      ICER = spinbscICER,
                      Labels = c("base case", Labels),
                      TypeRange = TypeRange,
                      FullName = c("Base case", rep(c(dsa.fullnames), each = dsa.levels)))

zolgspin <- data.frame(Scenario = 1:n.dsasims,
                       Parameter = c("base case", rep(dsa.rownames, each = dsa.levels)),
                       ScenarioNumb = c(1, rep(2:n, n.dsa)),
                       Bound = c("base case", rep(c(rep("lower", (dsa.levels)/2), rep("upper", dsa.levels/2)), n.dsa)),
                       ICER = zolgspinICER,
                       Labels = c("base case", Labels),
                       TypeRange = TypeRange,
                       FullName = c("Base case", rep(c(dsa.fullnames), each = dsa.levels)))

zolgbsc <- zolgbsc[zolgbsc$Parameter != "c.spin", ]
spinbsc <- spinbsc[spinbsc$Parameter != "c.zolg", ] 


# Save all generated data
write.csv2(zolgout, file.path("Output", "zolgoutDSA.csv"))
write.csv2(spinout, file.path("Output", "spinoutDSA.csv"))
write.csv2(bscout, file.path("Output", "bscoutDSA.csv"))
write.csv2(zolgbsc, file.path("Output", "zolgbscDSA.csv"))
write.csv2(spinbsc, file.path("Output", "spinbscDSA.csv"))
write.csv2(zolgspin, file.path("Output", "zolgspinDSA.csv"))

BRRR::skrrrahh(36)

