sDSA <- function(dat, MCSE, title, top = FALSE){
  # Set ICER values as numerical
  dat$ICER <- as.numeric(as.character(dat$ICER))
  # You may get a warning that NAs are introduced, this happens if ICERs are 
  # missing (e.g. not calculable ICERs)
  
  # Set NA ICERs at 0 (this default is decribed in the publication)
  dat$ICER[is.na(dat$ICER)] <- 0
  
  # Set the base case ICER
  ICER.basecase <- unique(dat$ICER[dat$Bound == 'base case'])
  
  # Calculate ICERs relative to the base case instead of absolute
  dat$relICER <- dat$ICER - ICER.basecase
  
  # Replace base case bound with NA value (necessary for plot generation)
  dat$Bound[dat$Bound == "base case"] <- NA
  
  # Create the names to be used in the plot, 
  # including the parameter name + the description of the range
  # dat$Parameter2 <- paste(dat$FullName, "(", dat$TypeRange, ")")
  # dat$Parameter2 <- gsub('\\(', '\n(', dat$Parameter2)
  dat$Parameter2 <- dat$FullName
  
  ## Construction of the order and the labels of the plot
  # Get a vector of the short parameter names
  parameters <- as.character(unlist(unique(dat$Parameter)))
  
  # Calculate the absolute range of each parameter (for plot order)
  absrange <- as.numeric(c(rep(NA, length(parameters))))
  for (i in 1:length(parameters)) {
    absrange[i] <- abs(max(dat$ICER[dat$Parameter == parameters[i]]) - 
                         min(dat$ICER[dat$Parameter == parameters[i]]))
  }
  
  dat$CatOrder <- c(0, rep(absrange[2:length(absrange)], each = max(dat$ScenarioNumb - 1)))
  orderdat <- order(dat$CatOrder)
  
  

  ## End of construction of the order and the labels of the plot
  
  
  ### ----------------End of getting data ready for the plot----------------- ###
  ### ----------------------------------------------------------------------- ###
  
  
  ### ----------------------------------------------------------------------- ###
  ### --------------------------Generating the plot-------------------------- ###
  
  # Set limits of the graph. With your own data you will probably have to try 
  # a few times before getting it right
  # At first try, a good shot of nice minimum and maximum values are those 
  # relative ICER values
  
  # Get minimum
  min(dat$relICER)
  # Get maximum
  max(dat$relICER)
  # Set limits (round down from min and up from max)
  xaxismin <- -60000
  xaxismax <- 60000
  limits <- c(xaxismin, xaxismax)
  
  # Set where labels are placed, 
  # the minimum should be slightly higher than the min limit
  breaks <- c(-100000, 100000)
  
  # Set how often labels are placed
  breaks.by <- 50000
  
  # Set the values displayed in the labels 
  # (should match the breaks & breaks_by arguments)
  #xlabels <- seq(-200000, 250000, 50000)
  
  # @ Rick, waarom dan niet:
  xlabels <- seq(breaks[1], breaks[2], breaks.by)
  
  # get top x values
  if (is.numeric(top)){  
    dat <- tail(dat[order(dat$CatOrder), ], top*14)
  }
  # sDSA plot
  gg1 <- ggplot(data = dat, aes(x = Parameter2, y = relICER, fill = Bound, order = order(Parameter2))) +
    geom_bar(stat = "identity", position = "dodge2",
             aes(x = fct_reorder(dat$Parameter2, dat$CatOrder, .desc = FALSE))) +
    theme(axis.ticks.y = element_blank()) +
    theme(panel.grid = element_blank()) +
    theme(panel.grid.major.x = element_line(colour = "grey80")) +
    theme(panel.grid.minor.x = element_line(colour = "grey80")) +
    scale_y_continuous(name = "Difference from base case ICER", labels = scales::number) +
    scale_fill_manual(breaks = c("lower", "upper"), 
                      values = c("#F8766D", "#00BFC4"),
                      label = c("Lower", "Upper")) +
    scale_x_discrete(name = "Parameter", 
                     breaks = dat$Parameter2, 
                     labels = dat$Parameter2) +
    theme(legend.title = element_blank(), 
          legend.position = "bottom", 
          legend.text = element_text(size = 56)) +
    theme(axis.text.x = element_text(size = 50), 
          axis.text.y = element_text(size = 52), 
          axis.title.x = element_text(size = 56), 
          axis.title.y = element_text(size = 56)) +
    geom_vline(xintercept = seq(1, length(dat[, 1]), 1) + 0.5, color = "grey74") +
    geom_hline(yintercept = 0, color = "grey50") +
    geom_hline(yintercept = -MCSE, color = "grey50", linetype = 2)+
    geom_hline(yintercept = MCSE, color = "grey50", linetype = 2)+
    coord_flip()+
    theme_Publication()+
    ggtitle(title)
  
  return(gg1)
}

