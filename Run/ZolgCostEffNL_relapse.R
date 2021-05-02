packs <- c("foreach", "doSNOW", "devtools", "ggplot2", "igraph")
for(i in 1:length(packs)){
  if (!packs[i] %in% installed.packages()) {install.packages(packs[i])}
}
if (!"dampack" %in% installed.packages()) {devtools::install_github("feralaes/dampack")}

library(foreach)
library(doSNOW)
library(ggplot2)

# Setup multicore cluster
cores <- parallel::detectCores()
cl <- makeSOCKcluster(cores[1] - 1)
registerDoSNOW(cl)

##################################### Model input ##########################################

# Model input
n.i   <- 10000 # number of simulated individuals
n.t   <- 1200                               # time horizon, 1200 cycles
v.n   <- c("1","2","3", "0",                   # Health states in the model
           "Death")
n.s   <- length(v.n)                           # the number of health states

d.c   <- 0.04                                  # discounting of costs by 4%
d.e   <- 0.015                                 # discounting of QALYs by 1.5%
set.seed(29)

p32r    <- p21r <- c(0, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05) 
times   <- rep(seq(12, 1188, 12), each = length(p32r))
p32r    <- rep(p32r, length(times)/length(p32r))
p21r    <- rep(p21r, length(times)/length(p21r))
n.rel   <- length(times)

v.dwc <- c(1, rep(1 / (1 + d.c) ^ (0:((n.t/12) -1)), each = 12))    # calculate the cost discount weight based on the discount rate d.c 
v.dwe <- c(1, rep(1 / (1 + d.e) ^ (0:((n.t/12) -1)), each = 12))    # calculate the QALY discount weight based on the discount rate d.e


# Input survival
p.1D    <- source(file.path("scripts", "p_1D.R"))[[1]]
p.2D    <- source(file.path("scripts", "p_2D.R"))[[1]]
p.3D    <- source(file.path("scripts", "lifetable_NL.R"))[[1]]
p.0D    <- source(file.path("scripts", "p_1D.R"))[[1]]             # probability to die in 0
p.10    <- rep(0, n.t)
p.0Dbsc    <- source(file.path("scripts", "p_1D.R"))[[1]]*1.35                    # probability to die in 0
p.10bsc    <- source(file.path("scripts", "p_10.R"))[[1]]


# Input costs and utilities
c.1     <- 9625*1.0125*1.0196
c.2     <- 10197*1.0125*1.0196
c.3     <- 5680*1.0125*1.0196
c.0     <- 14725.7*1.0125*1.0196


u.1     <- 0.733
u.2     <- 0.752
u.3     <- 0.878
u.0     <- 0.733


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

MicroSimRel <- function(n.i, n.t, v.n, d.c, d.e, f,
                     treat = c("Zolgensma", "BSC", "Spinraza"), relapse = FALSE, relapse.time, p.32r, p.21r) {
  # Arguments:  
  # v.M_1:    vector of initial states for individuals
  # n.i:      number of individuals
  # n.t:      total number of cycles to run the model
  # v.n:      vector of health state names 
  # d.c:      discount rate for costs
  # d.e:      discount rate for health outcome (QALYs)
  # TR.out:   should the output include a Microsimulation trace? (default is TRUE). Is used for trace plots with simplot()
  # TS.out:   should the output include a matrix of transitions between states? (default is TRUE)
  # data.out: should the output include the matrices with states, costs and effects (default is TRUE) (adds  +1 GB data in this model, only run when needed)
  # treat:    which treatment setting should be used for the model? Choice between Zolgensma and bscraza
  # seed:     starting seed number for random number generator (default is 1) 
  # Makes use of:
  # Probs:   function for the estimation of transition probabilities
  # Costs:   function for the estimation of cost state values
  # Effs:    function for the estimation of state specific health outcomes (QALYs)
  
  if (treat == "Zolgensma") {v.M_1 <- c(rep("1", n.i*0.08), rep("2", n.i*0.59), rep("3", n.i*0.33)); c.Trt <- 0; Trt <- TRUE}
  if (treat == "BSC")       {v.M_1 <- c(rep("1", n.i)); p.10 <- p.10bsc; p.0D <- p.0Dbsc; c.Trt <- 0; Trt <- FALSE}
  if (treat == "Spinraza")  {v.M_1 <- c(rep("1", n.i*0.5641), rep("2", n.i*0.1859), rep("3", n.i*0), rep("0", n.i*0.25)); c.Trt <- 83300; Trt <- TRUE}
  

  
  # Create the matrix capturing the state name/costs/health outcomes for all individuals at each time point 
  m.M <- m.C <- m.E <-  matrix(nrow = n.i, ncol = n.t + 1, 
                               dimnames = list(paste("ind", 1:n.i, sep = " "), 
                                               paste("cycle", 0:n.t, sep = " ")))  
  
  m.M[, 1] <- v.M_1                     # indicate the initial health state   
  
  m.C[, 1] <- Costs(m.M[, 1], Trt, t = 1, prev = NA, c.Trt)  # estimate costs per individual for the initial health state
  if (treat == "Zolgensma") {
    m.C[, 1] <- m.C[, 1] + 2001336.11
  }  # add €2000173 for treatment costs of Zolgensma   
  if (treat == "Spinraza"){
    m.C[, 1] <- m.C[, 1] + c.Trt*4 
  }
  
  m.E[, 1] <- Effs(m.M[, 1])  # estimate QALYs per individual for the initial health state
  
  reltime <- numeric(n.i)
  
  
  for (t in 1:n.t) {
    cur <- data.frame(m.M[, t], 1:n.i)
    cur <- cur[order(cur[, 1]), ]
    suppressWarnings(m.p <- Probs(cur, t, relapse, relapse.time, reltime, p.21r, p.32r))           # calculate the transition probabilities at cycle t 
    m.p <- m.p[order(cur[, 2]), ]
    m.M[, t + 1] <- samplev(prob = m.p[, -ncol(m.p)], m = 1)  # sample the next health state and store that state in matrix m.M 
    m.C[, t + 1] <- Costs(m.M[, t + 1], Trt, t = t, prev = m.M[, t], c.Trt)   # estimate costs per individual during cycle t + 1 conditional on treatment
    m.E[, t + 1] <- Effs(m.M[, t + 1])   # estimate QALYs per individual during cycle t + 1 conditional on treatment
    cat('\r', paste(round(t/n.t * 100, 2), "% done    ", sep = " "))       # display the progress of the simulation
    
    reltime <- ifelse((((m.M[, t + 1] != m.M[, t]) & m.M[, t + 1] != "Death") & m.M[, t + 1] != "0" ), t, reltime)
    
  } # close the loop for the time points 
  
  
  tc <- m.C %*% v.dwc            # total (discounted) cost per individual
  te <- (m.E/12) %*% v.dwe       # total (discounted) QALYs per individual 
  
  tc_hat <- mean(tc)       # average (discounted) cost 
  te_hat <- mean(te)        # average (discounted) QALYs
  
  results <- list(tc_hat = tc_hat, te_hat = te_hat) # store the results from the simulation in a list  
  return(results)  # return the results
  
  
}  # end of the MicroSim function  


#### Probability function
# The Probs function that updates the transition probabilities of every cycle is shown below.

Probs <- function(M_it, t = t, relapse, relapse.time, reltime, p.21r, p.32r) { 
  # M_it:    health state occupied by individual i at cycle t (character variable)
  
  m.p.it <- matrix(NA, n.s, n.i)     # create vector of state transition probabilities
  rownames(m.p.it) <- v.n  # assign vector names
  
  # update the v.p with the appropriate probabilities
  if (relapse == FALSE | (relapse == TRUE && t < relapse.time)){
    m.p.it[,M_it[, 1] == "1"]     <- c(1 - (p.1D[t] + p.10[t]), 0, 0, p.10[t], p.1D[t])                 
    m.p.it[,M_it[, 1] == "2"]     <- c(0, (1 - p.2D[t]), 0, 0, p.2D[t])   
    m.p.it[,M_it[, 1] == "3"]     <- c(0, 0, 1 - p.3D[t], 0, p.3D[t])                           
    m.p.it[,M_it[, 1] == "0"]     <- c(0, 0, 0, 1 - p.0D[t], p.0D[t])
    m.p.it[,M_it[, 1] == "Death"] <- c(0, 0, 0, 0, 1)
  }
  
  # Figure out exact position and length references for 2 and 3
  else if (relapse == TRUE && t >= relapse.time){
    m.p.it[,M_it[, 1] == "1"]     <- matrix(c(1 - (p.1D[t - reltime[1:length(m.p.it[1,M_it[, 1] == "1"])]] + p.10[t - reltime[1:length(m.p.it[1,M_it[, 1] == "1"])]]), rep(0, length(m.p.it[1,M_it[, 1] == "1"])), rep(0, length(m.p.it[1,M_it[, 1] == "1"])), p.10[t - reltime[1:length(m.p.it[1,M_it[, 1] == "1"])]], p.1D[t - reltime[1:length(m.p.it[1,M_it[, 1] == "1"])]]), nrow = n.s, byrow = TRUE)
    m.p.it[,M_it[, 1] == "2"]     <- matrix(c(rep(p.21r, length(m.p.it[2,M_it[, 1] == "2"])), 1 - (p.2D[t - reltime[(length(m.p.it[1,M_it[, 1] == "1"]) + 1):(length(m.p.it[1,M_it[, 1] == "1"]) + length(m.p.it[2,M_it[, 1] == "2"]))]] + p.21r), rep(0, length(m.p.it[2,M_it[, 1] == "2"])), rep(0, length(m.p.it[2,M_it[, 1] == "2"])), p.2D[t - reltime[(length(m.p.it[1,M_it[, 1] == "1"]) + 1):(length(m.p.it[1,M_it[, 1] == "1"]) + length(m.p.it[2,M_it[, 1] == "2"]))]]), nrow = n.s, byrow = TRUE)
    m.p.it[,M_it[, 1] == "3"]     <- matrix(c(rep(0, length(m.p.it[3,M_it[, 1] == "3"])), rep(p.32r, length(m.p.it[3,M_it[, 1] == "3"])), 1 - (p.3D[t - reltime[(length(m.p.it[2,M_it[, 1] == "2"]) + 1):(length(m.p.it[2,M_it[, 1] == "2"]) + length(m.p.it[3,M_it[, 1] == "3"]))]] + p.32r), rep(0, length(m.p.it[3,M_it[, 1] == "3"])), p.3D[t - reltime[(length(m.p.it[2,M_it[, 1] == "2"]) + 1):(length(m.p.it[2,M_it[, 1] == "2"]) + length(m.p.it[3,M_it[, 1] == "3"]))]]), nrow = n.s, byrow = TRUE)
    ifelse(length(m.p.it[,M_it[, 1] == "0"]) > 0, m.p.it[,M_it[, 1] == "0"]  <- matrix(c(rep(0, length(m.p.it[4,M_it[, 1] == "0"])), rep(0, length(m.p.it[4,M_it[, 1] == "0"])), rep(0, length(m.p.it[4,M_it[, 1] == "0"])), 1 - p.0D[t - reltime[(length(m.p.it[3,M_it[, 1] == "3"]) + 1):(length(m.p.it[3,M_it[, 1] == "3"]) + length(m.p.it[4,M_it[, 1] == "0"]))]], p.0D[t - reltime[(length(m.p.it[3,M_it[, 1] == "3"]) + 1):(length(m.p.it[3,M_it[, 1] == "3"]) + length(m.p.it[4,M_it[, 1] == "0"]))]]), nrow = n.s, byrow = TRUE),
                   m.p.it[,M_it[, 1] == "0"] <- c(0, 0, 0, 1 - p.0D[t], p.0D[t]))
    m.p.it[,M_it[, 1] == "Death"] <- c(0, 0, 0, 0, 1)
  }
  
  
  ifelse(colSums(m.p.it) == 1, return(matrix(c(t(m.p.it), M_it[, 2]), ncol = n.s + 1, dimnames = list(NULL, c(v.n, "ID")))), 
         print("Probabilities do not sum to 1"))   # return the transition probabilities or produce an error
  
  
}


### Costs function
# The Costs function estimates the costs at every cycle.

Costs <- function (M_it, Trt, t, prev, c.Trt) {
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
zolgout <- matrix(NA, n.rel, 2)
zolgtemp <- numeric(2)
bscout  <- matrix(NA, n.rel, 2)
bsctemp <- numeric(2)

# Setup progress bar
pb <- txtProgressBar(min = 1, max = n.rel, style = 3)
progress <- function(n) {setTxtProgressBar(pb, n)}
opts <- list(progress = progress)

p <- Sys.time()
p

zolgout <- foreach(f = 1:n.rel, .combine = rbind, .options.snow = opts) %dopar% {
  sim_zolg  <- MicroSimRel(n.i, n.t, v.n, d.c, d.e, treat = "Zolgensma",
                           f = f, relapse = TRUE, relapse.time = times[f], p.32r = p32r[f], p.21r = p21r[f]) # run for Zolgensma
  
  zolgtemp[1] <- sim_zolg$tc_hat
  zolgtemp[2] <- sim_zolg$te_hat
  
  setTxtProgressBar(pb, f)
  zolgtemp
}
close(pb)

pb <- txtProgressBar(min = 1, max = n.rel, style = 3)
bscout <- foreach(f = 1:n.rel, .combine = rbind, .options.snow = opts) %dopar% {
  sim_bsc  <- MicroSimRel(n.i, n.t, v.n, d.c, d.e, treat = "BSC",
                           f = f, relapse = TRUE, relapse.time = times[f], p.32r = p32r[f], p.21r = p21r[f]) # run for BSC
  
  bsctemp[1]  <- sim_bsc$tc_hat
  bsctemp[2]  <- sim_bsc$te_hat
  
  setTxtProgressBar(pb, f)
  bsctemp
}
close(pb)
stopCluster(cl)


zolgbsccost <- c(zolgout[, 1] - bscout[, 1])
zolgbsceff  <- c(zolgout[, 2] - bscout[, 2])
zolgbscICER <- zolgbsccost / zolgbsceff
zolgbsc <- data.frame(Costs = zolgbsccost, Effects = zolgbsceff, ICERs = zolgbscICER, relapseTime = times/12, relapseProb = p32r,
                      row.names = NULL)
comp.time <- Sys.time() - p
comp.time

zolgbsc

source(file.path("scripts", "gg_pubthemes.R"))
ggplot(zolgbsc, aes(relapseTime, ICERs, col = factor(relapseProb)))+
  geom_point()+
  scale_color_discrete(name = "Relapse prob.")+
  xlab("Year of relapse start")+
  theme_Publication()

write.csv2(zolgbsc, file.path("Output", "relapseICERs.csv"))
write.csv2(p32r, file.path("Output", "relapseProbs.csv"))
write.csv2(times, file.path("Output", "relapseTimes.csv"))
