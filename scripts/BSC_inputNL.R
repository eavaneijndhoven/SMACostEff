# Survival inputs
p.1D    <- source(here("scripts", "P_1D_BSC.R"))[[1]]         
p.2D    <- source(here("scripts", "p_2D.R"))[[1]]
p.3D    <- source(here("scripts", "lifetable_NL.R"))[[1]]
p.0D    <- source(here("scripts", "p_1D.R"))[[1]]                  
p.10    <- source(here("scripts", "p_10.R"))[[1]]


# Cost and utility inputs 
c.Trt   <- 0               # cost of treatment (per cycle)

c.1     <- 10701*1.0125*1.0196
c.2     <- 10457*1.0125*1.0196
c.3     <- 5808.5*1.0125*1.0196
c.0     <- 14725.7*1.0125*1.0196


u.1     <- 0.733
u.2     <- 0.752
u.3     <- 0.878
u.0     <- 0.733

## ICER UTILITIES
# u.1     <- 0.19
# u.2     <- 0.60
# u.3     <- c(rep(0.922, 360), rep(c(0.901, 0.871, 0.842, 0.823, 0.736), each = 120), rep(0.736, 240))
# u.0     <- 0.19

Trt     <- FALSE

