# Install packages (if needed) and load them
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

# Model input
set.seed(29)
n.i   <- 100                                # number of simulated individuals
n.t   <- 1200                               # time horizon, 30 cycles
v.n   <- c("1","2","3","0","Death")         # health state names
n.s   <- length(v.n)                        # the number of health states 
d.c   <- 0.04                               # discount rate costs
d.e   <- 0.015                              # discount rate effects

v.dwc <- c(1, rep(1 / (1 + d.c) ^ (0:((n.t/12) -1)), each = 12))    # calculate the cost discount weight based on the discount rate d.c 
v.dwe <- c(1, rep(1 / (1 + d.e) ^ (0:((n.t/12) -1)), each = 12))    # calculate the QALY discount weight based on the discount rate d.e

n.psa <- 1250                                  # number of PSA iterations

# Input treatment effect
source(file.path("scripts", "TreatPSA.R"))
TreatPSAzolg <- StatePSA(n.psa, "Zolgensma")
TreatPSAspin <- StatePSA(n.psa, "Spinraza")

# Input costs and utilities
source(file.path("scripts", "CostUtilPSA.R"))
costTrt   <- CostPSA(n.psa, Trt = TRUE)
costnoTrt <- CostPSA(n.psa, Trt = FALSE)
util      <- UtilPSA(n.psa)
zolg      <- TrtCost(n.psa, "Zolgensma")
spin      <- TrtCost(n.psa, "Spinraza")

# Input survival
source(file.path("scripts", "MicroPSAsurv.R"))   
p.1Dbsc <- p.1D.BSC.PSAfun(n.psa)
p.1D    <- p.1D.PSAfun(n.psa)
p.2D    <- p.2D.PSAfun(n.psa)
p.3D    <- p.3D.PSAfun(n.psa)
p.0D    <- p.1D
p.10bsc <- p.10.BSC.PSAfun(n.psa, p.1Dbsc)
p.10    <- matrix(0, n.t, n.psa)

p.2D    <- ifelse(p.2D > p.3D, p.2D, p.3D)


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
# The MicroSim function for the simple microsimulation of the 'Sick-Sicker' model keeps track of what happens to each individual during each cycle. 
MicroSimPSA <- function(n.i, n.t, v.n, d.c, d.e, 
                        treat = c("Zolgensma", "BSC"), f, zolg, cost) {
  
  if (treat == "Zolgensma") {v.M_1 <- c(rep("1", n.i*TreatPSAzolg[1, f] + 1), rep("2", n.i*TreatPSAzolg[2, f] + 1), rep("3", n.i*TreatPSAzolg[3, f] + 1)); Trt <- TRUE; Trtcost <- rep(0, n.psa)}
  if (treat == "BSC")       {v.M_1 <- c(rep("1", n.i*0.68), rep("0", n.i*0.32)); Trt <- FALSE; Trtcost <- rep(0, n.psa)}
  if (treat == "Spinraza")  {v.M_1 <- c(rep("1", n.i*TreatPSAspin[1, f] + 1), rep("2", n.i*TreatPSAspin[2, f] + 1), rep("0", n.i*TreatPSAspin[3, f] + 1)); Trt <- TRUE; Trtcost <- spin}
  
  
  # Create the matrix capturing the state name/costs/health outcomes for all individuals at each time point 
  m.M <- m.C <- m.E <-  matrix(nrow = n.i, ncol = n.t + 1, 
                               dimnames = list(paste("ind", 1:n.i, sep = " "), 
                                               paste("cycle", 0:n.t, sep = " ")))  
  
  m.M[, 1] <- v.M_1[1:n.i]                     # indicate the initial health state   
  
  
  m.C[, 1] <- Costs(m.M[, 1], f, prev = NA, t = 1, Trt = Trt, Trtcost, cost)  # estimate costs per individual for the initial health state
  m.C[, 1] <- m.C[, 1] + zolg[f] * Trt
  
  m.E[, 1] <- Effs (m.M[, 1], f)  # estimate QALYs per individual for the initial health state
  
  for (t in 1:n.t) {
    m.p <- Probs(m.M[, t], t = t, f = f, Trt)           # calculate the transition probabilities at cycle t 
    
    
    m.M[, t + 1] <- samplev(prob = m.p, m = 1)  # sample the next health state and store that state in matrix m.M 
    m.C[, t + 1] <- Costs(m.M[, t + 1], f, prev = m.M[, t], t = t, Trt = Trt, Trtcost, cost)   # estimate costs per individual during cycle t + 1 conditional on treatment
    m.E[, t + 1] <- Effs(m.M[, t + 1], f)   # estimate QALYs per individual during cycle t + 1 conditional on treatment
    
    
  } # close the loop for the time points 
  
  
  tc <- m.C %*% v.dwc       # total (discounted) cost per individual
  te <- (m.E/12) %*% v.dwe       # total (discounted) QALYs per individual 
  
  tc_hat <- mean(tc)        # average (discounted) cost 
  te_hat <- mean(te)        # average (discounted) QALYs
  
  results <- list(tc_hat = tc_hat, te_hat = te_hat) # store the results from the simulation in a list  
  return(results)  # return the results
  
  
}  # end of the MicroSim function  


#### Probability function
# The Probs function that updates the transition probabilities of every cycle is shown below.

Probs <- function(M_it, t, f, Trt) { 
  # M_it:    health state occupied by individual i at cycle t (character variable)
  
  m.p.it <- matrix(NA, n.s, n.i)     # create vector of state transition probabilities
  rownames(m.p.it) <- v.n            # assign names to the vector
  
  
  # update the v.p with the appropriate probabilities   
  ifelse(Trt == TRUE, m.p.it[,M_it == "1"]     <- c(1 - (p.1D[t, f] + p.10[t, f]), 0, 0, p.10[t, f], p.1D[t, f]),
         m.p.it[,M_it == "1"]     <- c(1 - (p.1Dbsc[t, f] + p.10bsc[t, f]), 0, 0, p.10bsc[t, f], p.1Dbsc[t, f]))
  m.p.it[,M_it == "2"]     <- c(0, (1 - p.2D[t, f]), 0, 0, p.2D[t, f])   
  m.p.it[,M_it == "3"]     <- c(0, 0, 1 - p.3D[t, f], 0, p.3D[t, f])                           
  m.p.it[,M_it == "0"]     <- c(0, 0, 0, 1 - p.0D[t, f], p.0D[t, f])
  m.p.it[,M_it == "Death"] <- c(0, 0, 0, 0, 1)
  
  ifelse(colSums(m.p.it) == 1, return(t(m.p.it)), print("Probabilities do not sum to 1"))   # return the transition probabilities or produce an error
  
  
  
}


### Costs function
# The Costs function estimates the costs at every cycle.

Costs <- function (M_it, f, prev, t, Trt, Trtcost, cost) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  # Trt:  is the individual being treated? (default is FALSE) 
  
  c.it <- 0                                  # by default the cost for everyone is 0
  c.it[M_it == "1"]     <- ifelse(t %in% seq(5, 1200, by = 4),  cost[[1]][f] + Trtcost[f] * Trt, cost[[1]][f])               
  c.it[M_it == "2"]     <- ifelse(t %in% seq(5, 1200, by = 4),  cost[[2]][f] + Trtcost[f] * Trt, cost[[2]][f])
  c.it[M_it == "3"]     <- ifelse(t %in% seq(5, 1200, by = 4),  cost[[3]][f] + Trtcost[f] * Trt, cost[[3]][f])
  c.it[M_it == "0"]     <- ifelse(t %in% seq(5, 1200, by = 4),  cost[[4]][f] + Trtcost[f] * Trt, cost[[4]][f])               
  c.it[M_it == "Death"] <- 0
  c.it[M_it == "Death" & prev == "1"] <- cost[[4]][f]*3    # dying gives a one time cost      		                   # return the costs
  
  return(c.it)
}


### Health outcome function 
# The Effs function to update the utilities at every cycle.

Effs <- function (M_it, f, cl = 1) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  # Trt:  is the individual treated? (default is FALSE) 
  # cl:   cycle length (default is 1)
  
  u.it <- 0                      # by default the utility for everyone is zero
  u.it[M_it == "1"]     <- util[[1]][f]      
  u.it[M_it == "2"]     <- util[[2]][f] 
  u.it[M_it == "3"]     <- util[[3]][f]
  u.it[M_it == "0"]     <- util[[4]][f]      
  u.it[M_it == "Death"] <- 0
  
  QALYs <-  u.it * cl            # calculate the QALYs during cycle t
  return(QALYs) 
  
} 

# Setup progress bar
pb <- txtProgressBar(min = 1, max = n.psa, style = 3)
progress <- function(n) {setTxtProgressBar(pb, n)}
opts <- list(progress = progress)

# Setup output files
zolgout <- matrix(NA, n.psa, 2)
zolgtemp <- numeric(2)
spinout <- matrix(NA, n.psa, 2)
spintemp <- numeric(2)
bscout  <- matrix(NA, n.psa, 2)
bsctemp <- numeric(2)


# Run PSA for each treatment type separately. Multicore approach using foreach() to speed up process.
p <- Sys.time()
p
zolgout <- foreach(f = 1:n.psa, .combine = rbind, .options.snow = opts) %dopar% {
    sim_zolg  <- MicroSimPSA(n.i, n.t, v.n, d.c, d.e, treat = "Zolgensma",
                             f = f, zolg = zolg, cost = costTrt) # run for Zolgensma
    
    zolgtemp[1] <- sim_zolg$tc_hat
    zolgtemp[2] <- sim_zolg$te_hat
    
    setTxtProgressBar(pb, f)
    zolgtemp
  }


pb <- txtProgressBar(min = 1, max = n.psa, style = 3)
bscout <- foreach(f = 1:n.psa, .combine = rbind, .options.snow = opts) %dopar% {
    sim_bsc   <- MicroSimPSA(n.i, n.t, v.n, d.c, d.e, treat = "BSC",
                             f = f, zolg = rep(0, n.psa), cost = costnoTrt)  # run for BSC
    bsctemp[1]  <- sim_bsc$tc_hat
    bsctemp[2]  <- sim_bsc$te_hat
    
    setTxtProgressBar(pb, f)
    bsctemp
  }
close(pb)

pb <- txtProgressBar(min = 1, max = n.psa, style = 3)
spinout <- foreach(f = 1:n.psa, .combine = rbind, .options.snow = opts) %dopar% {
    sim_spin  <- MicroSimPSA(n.i, n.t, v.n, d.c, d.e, treat = "Spinraza",
                             f = f, zolg = rep(0, n.psa), cost = costTrt)
    spintemp[1] <- sim_spin$tc_hat
    spintemp[2] <- sim_spin$te_hat
    
    setTxtProgressBar(pb, f)
    spintemp
  }
close(pb)
stopCluster(cl)


# Combine results of PSAs to calculate incremental costs and effects and ICERs
zolgbsccost <- c(zolgout[, 1] - bscout[, 1])
zolgbsceff  <- c(zolgout[, 2] - bscout[, 2])
zolgbscICER <- zolgbsccost / zolgbsceff
zolgbsc <- matrix(c(zolgbsccost, zolgbsceff, zolgbscICER),
                  nrow = n.psa, ncol = 3)

spinbsccost <- c(spinout[, 1] - bscout[, 1])
spinbsceff  <- c(spinout[, 2] - bscout[, 2])
spinbscICER <- spinbsccost / spinbsceff
spinbsc <- matrix(c(spinbsccost, spinbsceff, spinbscICER),
                  nrow = n.psa, ncol = 3)


zolgspincost <- c(zolgout[, 1] - spinout[, 1])
zolgspineff  <- c(zolgout[, 2] - spinout[, 2])
zolgspinICER <- zolgspincost / zolgspineff
zolgspin <- matrix(c(zolgspincost, zolgspineff, zolgspinICER),
                   nrow = n.psa, ncol = 3)


comp.time <- Sys.time() - p
comp.time

####### RUN THIS AT FINAL MODEL


# Save all generated data
# write.csv2(p.1D, file.path("Output", "PSA", "p_1D_prob.csv"))
# write.csv2(p.2D, file.path("Output", "PSA", "p_2D_prob.csv"))
# write.csv2(p.3D, file.path("Output", "PSA", "p_3D_prob.csv"))
# write.csv2(p.0D, file.path("Output", "PSA", "P_0D_prob.csv"))
# write.csv2(p.1Dbsc, file.path("Output", "PSA", "P_1Dbsc_prob.csv"))
# write.csv2(p.10bsc, file.path("Output", "PSA", "P_10_prob.csv"))
# write.csv2(as.data.frame(costTrt), file.path("Output", "PSA", "costTrt_prob.csv"))
# write.csv2(as.data.frame(costnoTrt), file.path("Output", "PSA", "costnoTrt_prob.csv"))
# write.csv2(as.data.frame(util), file.path("Output", "PSA", "util_prob.csv"))
# write.csv2(zolg, file.path("Output", "PSA", "zolg_cost_prob.csv"))
# write.csv2(spin, file.path("Output", "PSA", "spin_cost_prob.csv"))
# write.csv2(TreatPSAzolg, file.path("Output", "PSA", "zolg_eff.csv"))
# write.csv2(TreatPSAspin, file.path("Output", "PSA", "spin_eff.csv"))
# write.csv2(zolgbsc, file.path("Output", "PSA", "zolgbsc.csv"))
# write.csv2(spinbsc, file.path("Output", "PSA", "spinbsc.csv"))
# write.csv2(zolgspin, file.path("Output", "PSA", "zolgspin.csv"))
# write.csv2(zolgout, file.path("Output", "PSA", "zolgout.csv"))
# write.csv2(spinout, file.path("Output", "PSA", "spinout.csv"))
# write.csv2(bscout, file.path("Output", "PSA", "bscout.csv"))


# Generate dataset for CE plane
source(file.path("scripts", "gg_pubthemes.R"))
zolgbscCE  <- cbind(zolgbsc[, 1:2], set = rep(1, nrow(zolgbsc)), shp = rep(1, nrow(zolgbsc)))
colnames(zolgbscCE) <- c("costs", "effects", "set", "shp")
zolgbscCE <- rbind(zolgbscCE, data.frame(costs = c(mean(zolgbsc[, 1]), 3106038), 
                                         effects = c(mean(zolgbsc[, 2]), 22.381),
                                         set = c(4, 7),
                                         shp = c(2, 2)))

spinbscCE  <- cbind(spinbsc[, 1:2], set = rep(2, nrow(spinbsc)), shp = rep(1, nrow(spinbsc)))
colnames(spinbscCE) <- c("costs", "effects", "set", "shp")
spinbscCE  <- rbind(spinbscCE, data.frame(costs = c(mean(spinbsc[, 1]), 2082665), 
                                          effects = c(mean(spinbsc[, 2]), 3.236),
                                          set = c(5, 8),
                                          shp = c(2, 2)))

zolgspinCE <- cbind(zolgspin[, 1:2], set = rep(3, nrow(zolgspin)), shp = rep(1, nrow(zolgspin)))
colnames(zolgspinCE) <- c("costs", "effects", "set", "shp")
zolgspinCE <- rbind(zolgspinCE, data.frame(costs = c(mean(zolgspin[, 1]), 1025996), 
                                           effects = c(mean(zolgspin[, 2]), 19.152),
                                           set = c(6, 9),
                                           shp = c(2, 2)))

ceplane <- rbind(zolgbscCE, spinbscCE, zolgspinCE)
colnames(ceplane) <- c("Costs", "Effects", "set", "shp")

# Make and save the CE plane
ceplaneplotfull <- ggplot(ceplane, aes(Effects, Costs, col = factor(set), shape = factor(shp), size = factor(shp), alpha = factor(shp)))+
  xlim(-30, 30)+
  xlab("Effects (QALYs)")+
  scale_y_continuous(name = "Costs (€)",
                     breaks = seq(from = -6000000, to = 6000000, by = 2000000),
                     labels = format(seq(from = -6000000, to = 6000000, by = 2000000), scientific = FALSE),
                     limits = c(-6000000, 6000000))+
  geom_abline(aes(intercept = 0, slope = 80000, color = "green"), lty = 5, show.legend = TRUE)+
  geom_abline(aes(intercept = 0, slope = 150000, color = "purple"), lty = 5, show.legend = TRUE)+
  geom_abline(aes(intercept = 0, slope = 300000, color = "slateblue"), lty = 5, show.legend = TRUE)+
  geom_hline(yintercept = 0, color = "gray25", lty = 3)+
  geom_vline(xintercept = 0, color = "gray25", lty = 3)+
  geom_point()+
  scale_colour_manual(name = "Legend",
                      labels = c("Prob. ICER Zolg/BSC", "Prob. ICER Spin/BSC", "Prob. ICER Zolg/Spin",
                                 "Prob. mean Zolg/BSC", "Prob. mean Spin/BSC", "Prob. mean Zolg/Spin",
                                 "Det. mean Zolg/BSC", "Det. mean Zolg/BSC", "Det. mean Zolg/Spin",
                                 "€80000/QALY", "€150000/QALY", "€300000/QALY"),
                      values = c("black", "red", "orange",
                                 "gray70", "pink", "lightsalmon3",
                                 "lightblue", "firebrick4", "darkorange4",
                                 "green", "purple", "slateblue"),
                      guide = guide_legend(override.aes = list(
                        linetype = c(rep(0, 9), rep(5, 3)),
                        shape = c(rep(16, 3), rep(17, 6), rep(NA, 3))
                      )))+
  scale_shape_discrete(guide = FALSE)+
  scale_size_discrete(guide = FALSE, range = c(0.9, 3))+
  scale_alpha_discrete(guide = FALSE, range = c(0.5, 1))+
  ggtitle("CE Plane of Zolgensma vs. Spinraza vs. BSC")+
  theme_Publication()

# ggsave(file.path("Output", "PSA", "CEplane.tiff"), plot = ceplaneplotfull, device = "tiff",
#        width = 20, height = 20, units = "cm", dpi = 320)

ceplaneplotfull


BRRR::skrrrahh(36)