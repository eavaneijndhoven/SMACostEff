# Survival inputs
p.1D    <- source(here("scripts", "p_1D.R"))[[1]]
p.2D    <- source(here("scripts", "p_2D.R"))[[1]]
p.3D    <- source(here("scripts", "lifetable_NL.R"))[[1]]
p.0D    <- source(here("scripts", "p_1D.R"))[[1]]             # probability to die in 0
p.10    <- rep(0, n.t)
p.2D    <- ifelse(p.2D > p.3D, p.2D, p.3D)

# Cost and utility inputs 
c.Trt   <- 0               # cost of treatment (per cycle)

c.1     <- 9625*1.0125*1.0196
c.2     <- 10197*1.0125*1.0196
c.3     <- 5680*1.0125*1.0196
c.0     <- 14725.7*1.0125*1.0196


u.1     <- 0.733
u.2     <- 0.752
u.3     <- 0.878
u.0     <- 0.733


# ICER UTILITIES
# u.1     <- 0.29
# u.2     <- 0.65
# u.0     <- 0.19
# u.3     <- c(rep(0.922, 360), rep(c(0.901, 0.871, 0.842, 0.823, 0.736), each = 120), rep(0.736, 240))

Trt     <- TRUE




