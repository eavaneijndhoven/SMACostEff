simplot <- function(sim1, sim2, type = "area"){
  ## ggplot2 is used to make the plot
  library(ggplot2)
  
  # Extract trace plots and treatment names from MicroSim output
  sims <- list(as.data.frame(sim1$TR), as.data.frame(sim2$TR))
  simnames <- c(sim1$treat, sim2$treat)
  
  # Convert TR dataframes to long format
  sims <- lapply(sims, function(x){
    x$cycle <- 0:n.t
    x <- reshape(x, idvar = "cycle", varying = v.n,
                 v.names = "state", direction = "long")
    x$time <- factor(rep(v.n, each = n.t + 1))
    x$time <- factor(x$time, levels = rev(levels(x$time)))
    return(x)
  })
  
  # Add column with treatment name to TR, and rbind together
  sims <- mapply(function(x, y) cbind(x, treat = y), sims, c(simnames), SIMPLIFY = FALSE)  
  simplot <- do.call(rbind, sims)
  
  
  # Make the plots as stacked area plot, or lines
  if (type == "area"){  
    x <- ggplot(data = simplot, aes(cycle, state, fill = time, order = as.numeric(time)))+
      geom_area()+
      ylab("Proportion")+
      scale_fill_discrete(name = "Health state")+
      facet_wrap(~ treat, nrow = 2)+
      theme_classic()
  }
  
  if (type == "line"){
    x <- ggplot(data = simplot, aes(cycle, state, col = time, order = as.numeric(time)))+
      geom_line(size = 1.5)+
      ylab("Proportion")+
      scale_color_discrete(name = "Health state")+
      facet_wrap(~ treat, nrow = 2)+
      theme_classic()
  }
  
  # Make the plot
  return(x)
}

