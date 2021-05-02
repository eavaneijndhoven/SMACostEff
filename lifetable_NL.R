# Read data
cbs <- read.csv2(here("Files", "Overlevingskansen__geslacht__leeftijd_28102019_115617.csv"))

# Convert percentages into proportions
cbs$Survival <- cbs$Survival/100

# Get mortality from survival
cbs$Mortality <- 1 - cbs$Survival

# Calculate montly probabilities from yearly mortality
monthrate <- -1/12*log(1-cbs$Mortality)
monthtp <- 1-exp(-monthrate)
monthtp <- rep(monthtp, each = 12)

# Repeat survival 98-99 for 99-100 interval
monthtp[1189:1200] <- monthtp[1177:1188]
monthtp
