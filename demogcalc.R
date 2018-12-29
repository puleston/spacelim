# function demogcalc
# calculates demographic stats from survival and fertility vectors

demogcalc <- function(px,mx){
  
  # probability of dying in each interval
  qx <- (1-px)

  ################################### calc life expectancy
  
  # calculate cumulative survivorship (prob of living to age x)
  lx <- rep(0,length(px))
  lx[1] <- 1
  for(i in 2:length(lx)){
    lx[i] <- lx[(i-1)]*px[(i-1)]
  }
  
  # ax = person-years lived by people who died in interval
  # 1st entry captures infant mort, fixed, unless mort is low; females, adapted from Preston et al., p. 48
  ax <- rep(1,length(px))/2
  ax[1] <- min(0.35,(0.053 + 2.8*qx[1]));
  
  # Person-years lived in the interval by survivors and those who died
  Lx <- rep(0,length(px))
  for(i in 1:length(Lx)){
  Lx[i] <- lx[i+1] + ax[i]*(lx[i]-lx[i+1])
  }
  
# person-years remaining to those just born (life expectancy at birth)  
e0 <- sum(Lx[1:(length(Lx)-1)])

# remaining life expectancy of those who make it to age 10
e10 <- sum(Lx[11:(length(Lx)-1)])/lx[11]

################################### calc fertility rates

# Total fertility rate (here, number of daughters born to a baby girl who lives through reproductive years)
TFR <- sum(mx)

# Net reproductive rate (number of daughters a female baby is expected to have after accounting for her mortality)
NRR <- sum((lx*mx))

################################### calc MAC or generation time

xx <- c(1:length(px))

muC <- sum((lx*mx*(xx-0.5)))/NRR

################################### return results

list(e0,e10,TFR,NRR,muC)
}


