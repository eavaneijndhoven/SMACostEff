# Randomize health states for treatment effect using dirichlet distribution based on base case health states
StatePSA <- function(n.psa, Treat = c("Zolgensma", "Spinraza")){
  
  ifelse(Treat == "Zolgensma", psa.state <- list(state.1 = 0.08, state.2 = 0.59, state.3 = 0.33),
                               psa.state <- list(state.1 = 0.5641, state.2 = 0.1849, state.0 = 0.25))
  
  psa.state <- lapply(psa.state, function (x) {
    append(x, (x*1.5-x*0.5)/(2*1.96))
    
  })
  
  psa.c <- dampack::dirichlet_params(p.mean = unlist(psa.state)[c(1, 3, 5)], sigma = unlist(psa.state)[c(2, 4, 6)])
  psa.out <- igraph::sample_dirichlet(n.psa, psa.c)
  
  return(psa.out)
}
