a <- Sys.time()

# Here library is used to locate all files. 
# You can also use base R's file.path(), but this makes sure connections also work in Rmarkdown
library(here)
##################################### Model input ##########################################

# Model input
n.i   <- 10000                                 # number of simulated individuals
n.t   <- 1200                                  # time horizon, 1200 cycles
v.n   <- c("1","2","3", "0",                   # Health states in the model
           "Death")
n.s   <- length(v.n)                           # the number of health states

d.c   <- 0.04                                  # discounting of costs by 4%
d.e   <- 0.015                                 # discounting of QALYs by 1.5%
set.seed(29)                                   # set the seed for reproducible randomization


##################################### Functions ###########################################

# THE NEW samplev() FUNCTION
# efficient implementation of the rMultinom() function of the Hmisc package #### 
# this function calculates health states for the next cycle based on individual probabilities

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

MicroSim <- function(n.i, n.t, v.n, d.c, d.e, TR.out = TRUE, TS.out = TRUE, data.out = TRUE,
                     treat = c("Zolgensma", "BSC", "Spinraza")) {
  # Arguments:  
  # 
  # n.i:      number of individuals
  # n.t:      total number of cycles to run the model
  # v.n:      vector of health state names 
  # d.c:      discount rate for costs
  # d.e:      discount rate for health outcome (QALYs)
  # TR.out:   should the output include a Microsimulation trace? (default is TRUE). Is used for trace plots with simplot()
  # TS.out:   should the output include a matrix of transitions between states? (default is TRUE)
  # data.out: should the output include the matrices with states, costs and effects (default is TRUE) (adds ~1 GB data in this model, only run when needed)
  # treat:    which treatment setting should be used for the model? Choice between Zolgensma and Spinraza
  # seed:     starting seed number for random number generator (default is 1) 
  # Makes use of:
  # v.M_1:    vector of initial states for individuals
  # Probs:   function for the estimation of transition probabilities
  # Costs:   function for the estimation of cost state values
  # Effs:    function for the estimation of state specific health outcomes (QALYs)
   
  # Depending on the treatment setting simulated, load costs, utilities and survival probabilities. v.M_1 places individuals in health states
  # according to treatment effect
  if (treat == "Zolgensma") {source(here("scripts", "zolginput_NL.R")); v.M_1 <- c(rep("1", n.i*0.08), rep("2", n.i*0.59), rep("3", n.i*0.33))}
  if (treat == "BSC")       {source(here("scripts", "BSC_inputNL.R")); v.M_1 <- c(rep("1", n.i*0.68), rep("0", n.i*0.32))}
  if (treat == "Spinraza")  {source(here("scripts", "spininput_NL.R")); v.M_1 <- c(rep("1", n.i*0.5641), rep("2", n.i*0.1859), rep("3", n.i*0), rep("0", n.i*0.25))}
  
  v.dwc <- c(1, rep(1 / (1 + d.c) ^ (0:((n.t/12) -1)), each = 12))    # calculate the cost discount weight based on the discount rate d.c 
  v.dwe <- c(1, rep(1 / (1 + d.e) ^ (0:((n.t/12) -1)), each = 12))    # calculate the QALY discount weight based on the discount rate d.e
  
  # Create the matrix capturing the state name/costs/health outcomes for all individuals at each time point 
  m.M <- m.C <- m.E <-  matrix(nrow = n.i, ncol = n.t + 1, 
                               dimnames = list(paste("ind", 1:n.i, sep = " "), 
                                               paste("cycle", 0:n.t, sep = " ")))  
  
  m.M[, 1] <- v.M_1                     # indicate the initial health state   
  
  m.C[, 1] <- Costs(m.M[, 1], Trt, t = 1, prev = NA)  # estimate costs per individual for the initial health state
  
  # add €2001336.11 for treatment costs of Zolgensma   
  if (treat == "Zolgensma") {
    m.C[, 1] <- m.C[, 1] + 2001336.11
  } 
  # add costs for four loading doses of Spinraza
  if (treat == "Spinraza"){
    m.C[, 1] <- m.C[, 1] + c.Trt*4 
  }
  
  m.E[, 1] <- Effs(m.M[, 1])  # estimate QALYs per individual for the initial health state
  
  for (t in 1:n.t) {
    m.p <- Probs(m.M[, t], t)           # calculate the transition probabilities at cycle t 
    
    
    m.M[, t + 1] <- samplev(prob = m.p, m = 1)  # sample the next health state and store that state in matrix m.M 
    m.C[, t + 1] <- Costs(m.M[, t + 1], Trt, t = t, prev = m.M[, t])   # estimate costs per individual during cycle t + 1 conditional on treatment
    m.E[, t + 1] <- Effs(m.M[, t + 1])   # estimate QALYs per individual during cycle t + 1 conditional on treatment
    cat('\r', paste(round(t/n.t * 100, 2), "% done    ", sep = " "))       # display the progress of the simulation
    
    
  } # close the loop for the time points 
  
  
  tc <- m.C %*% v.dwc            # total (discounted) cost per individual
  te <- (m.E/12) %*% v.dwe       # total (discounted) QALYs per individual. Monthly cycles and yearly utilities in input, so divide by 12 

  tc_hat <- mean(tc)       # average (discounted) cost 
  te_hat <- mean(te)        # average (discounted) QALYs
  
  if (TS.out == TRUE) {  # create a matrix of transitions across states
    TS <- paste(m.M, cbind(m.M[, -1], NA), sep = "->") # transitions from one state to the other
    TS <- matrix(TS, nrow = n.i)
    rownames(TS) <- paste("Ind",   1:n.i, sep = " ")   # name the rows 
    colnames(TS) <- paste("Cycle", 0:n.t, sep = " ")   # name the columns 
  } else {
    TS <- NULL
  }
  
  if (TR.out == TRUE) {
    TR <- t(apply(m.M, 2, function(x) table(factor(x, levels = v.n, ordered = TRUE))))
    TR <- TR / n.i                                       # create a distribution trace
    rownames(TR) <- paste("Cycle", 0:n.t, sep = " ")     # name the rows 
    colnames(TR) <- v.n                                  # name the columns
  } else {
    TR <- NULL
  }
  if (data.out == TRUE) {results <-  list(m.M = m.M, m.C = m.C, m.E = m.E, tc = tc, te = te, tc_hat = tc_hat, te_hat = te_hat, TS = TS, TR = TR, treat = treat)}
  else {results <-  list(tc = tc, te = te, tc_hat = tc_hat, te_hat = te_hat, TS = TS, TR = TR, treat = treat)} # store the results from the simulation in a list  
  return(results)  # return the results
  
  
  
}  # end of the MicroSim function  


#### Probability function
# The Probs function that updates the transition probabilities of every cycle is shown below.

Probs <- function(M_it, t = t) { 
  # M_it:    health state occupied by individual i at cycle t (character variable)
  
  m.p.it <- matrix(NA, n.s, n.i)     # create vector of state transition probabilities
  rownames(m.p.it) <- v.n            # assign names to the vector
  
  # update the v.p with the appropriate probabilities
  m.p.it[,M_it == "1"]     <- c(1 - (p.1D[t] + p.10[t]), 0, 0, p.10[t], p.1D[t])                 
  m.p.it[,M_it == "2"]     <- c(0, (1 - p.2D[t]), 0, 0, p.2D[t])   
  m.p.it[,M_it == "3"]     <- c(0, 0, 1 - p.3D[t], 0, p.3D[t])                           
  m.p.it[,M_it == "0"]     <- c(0, 0, 0, 1 - p.0D[t], p.0D[t])
  m.p.it[,M_it == "Death"] <- c(0, 0, 0, 0, 1)
  
  
  ifelse(colSums(m.p.it) == 1, return(t(m.p.it)), print("Probabilities do not sum to 1"))   # return the transition probabilities or produce an error
  
  
  
}


### Costs function
# The Costs function estimates the costs at every cycle.

Costs <- function (M_it, Trt, t, prev) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  
  
  c.it <- 0                                  # by default the cost for everyone is 0
  c.it[M_it == "1"]     <- ifelse(t %in% seq(5, 1200, by = 4),  c.1 + c.Trt * Trt, c.1)               
  c.it[M_it == "2"]     <- ifelse(t %in% seq(5, 1200, by = 4),  c.2 + c.Trt * Trt, c.2)
  c.it[M_it == "3"]     <- ifelse(t %in% seq(5, 1200, by = 4),  c.3 + c.Trt * Trt, c.3)
  c.it[M_it == "0"]     <- c.0                
  c.it[M_it == "Death"] <- 0
  c.it[M_it == "Death" & prev == "1"] <- c.0*3    # dying gives a one time cost
  
  return(c.it)        		                   # return the costs
  
}


### Health outcome function 
# The Effs function to update the utilities at every cycle.

Effs <- function (M_it, cl = 1, t) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  # Trt:  is the individual treated? (default is FALSE) 
  # cl:   cycle length (default is 1)
  
  u.it <- 0                      # by default the utility for everyone is zero
  u.it[M_it == "1"]     <- u.1      
  u.it[M_it == "2"]     <- u.2 
  u.it[M_it == "3"]     <- u.3
  u.it[M_it == "0"]     <- u.0      
  u.it[M_it == "Death"] <- 0
  
  QALYs <-  u.it * cl            # calculate the QALYs during cycle t
  return(QALYs)                  # return the QALYs
  
}


##################################### Run the simulation ##################################
# START SIMULATION
# At treat: input which two treatment settings to compare.
p = Sys.time()
sim_1     <- MicroSim(n.i, n.t, v.n, d.c, d.e, treat = "BSC", TS.out = FALSE, TR.out = TRUE, data.out = TRUE) 
sim_2     <- MicroSim(n.i, n.t, v.n, d.c, d.e, treat = "Spinraza", TS.out = FALSE, TR.out = TRUE, data.out = TRUE)  
comp.time <- Sys.time() - p

comp.time

################################# Cost-effectiveness analysis #############################
# store the mean costs (and the MCSE) of each strategy in a new variable C (vector costs)
v.C  <- c(sim_1$tc_hat, sim_2$tc_hat) 
sd.C <- c(sd(sim_1$tc), sd(sim_2$tc)) / sqrt(n.i)
# store the mean QALYs (and the MCSE) of each strategy in a new variable E (vector effects)
v.E  <- c(sim_1$te_hat, sim_2$te_hat)
sd.E <- c(sd(sim_1$te), sd(sim_2$te)) / sqrt(n.i)

delta.C <- v.C[2] - v.C[1]                   # calculate incremental costs
delta.E <- v.E[2] - v.E[1]                   # calculate incremental QALYs
sd.delta.E <- sd(sim_2$te - sim_1$te) / sqrt(n.i) # Monte Carlo Squared Error (MCSE) of incremental costs
sd.delta.C <- sd(sim_2$tc - sim_1$tc) / sqrt(n.i) # Monte Carlo Squared Error (MCSE) of incremental QALYs
ICER    <- delta.C / delta.E                 # calculate the ICER
results <- c(delta.C, delta.E, ICER)         # store the values in a new variable



# Create full incremental cost-effectiveness analysis table
table_micro <- data.frame(
  c(round(v.C, 0),  ""),           # costs per arm
  c(round(sd.C, 0), ""),           # MCSE for costs
  c(round(v.E, 3),  ""),           # health outcomes per arm
  c(round(sd.E, 3), ""),           # MCSE for health outcomes
  c("", round(delta.C, 0),   ""),  # incremental costs
  c("", round(sd.delta.C, 0),""),  # MCSE for incremental costs
  c("", round(delta.E, 3),   ""),   # incremental QALYs 
  c("", round(sd.delta.E, 3),""),  # MCSE for health outcomes (QALYs) gained
  c("", round(ICER, 0),      "")   # ICER
)
rownames(table_micro) <- c(c(sim_1$treat, sim_2$treat), "* are MCSE values")  # name the rows
colnames(table_micro) <- c("Costs", "*",  "QALYs", "*", "Incremental Costs", "*", "QALYs Gained", "*", "ICER") # name the columns


# Check transitions
# write.csv2(sim_2$TS, file.path("Output", "relapse2.csv"))


# Make trace plot
source(file.path("scripts", "simplotv2.R"))
simplot(sim_1, sim_2)
table_micro  # print the table

Sys.time() - a