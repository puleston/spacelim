
###################################################################################################################################
######################################################################################################### Contact and Personal Info

#   Food Limited Demography Run File
#   Code based on models described in Lee & Tuljapurkar (2008), Puleston & Tuljapurkar (2008) and Lee et al. (2009)
#   R Interpetation by Cedric Puleston and Cody Ross (ctross at ucdavis dot edu)
#   Major Version 1
#   Minor Version 3b

###################################################################################################################################
############################################################################################################ Parameter Definitions

# Supplied parameters:
# p0 - vector - survival prob from age class to next at full food (only use when E>=1)
# m0 - vector - fertility prob from age class to next at full food (only use when E>=1)
# pstar - vector - baseline survival prob for use w E < 1, they make gamma cdf give correct answer when E=1
# mstar - vector - baseline fertility prob for use w E < 1, they make gamma cdf give correct answer when E=1
# mortparms - vector - alpha of gamma cdf, responsivness of survival on food
# fertparm - scalar - alpha of gamma cdf, responsiveness of fertility to food
# rho - vector - consumption for each age class, relative to highest consuming age class
# phi - vector - worker effectivness, relative to most effective age class
# N.AC - scalar, number of age classes
# Y - scalar, kcal/ha/day in ag yield
# Y.CV - scalar, coeff of variation of Y
# Yac - scalar, annual serial autocorrelation of Y
# Am - scalar, total area of ag land available for cultivation in ha
# h - scalar capturing conversion between labor availability and land in cultivation (wo considering interference competition)
# J - scalar representing the daily caloric need of the calorically neediest age class

# Machinery parameters
# N.guess - scalar to get fzero function started in the right neighborhood for Nhat
# span - scalar indicating the length of time series in years

# Derived parameters:
# phi1 - scalar representing equilibrium age structure-weighted per capita work effectiveness
# rho1 - scalar representing equilibrium age structure-weighted per capita consumption

###################################################################################################################################
############################################################################################################ Load Parameters

SS<- "C:/Users/Empyrean/Desktop/Dropbox/Food Limited Demography and the Evolution of Inequality"
# Set to your local directory, where the .csv file of model parameters is stored
setwd(SS)

# Read data and store as d
d <- read.csv("demog_vectors.csv")

#################################################################################################### Set Pop1 parms
# Assign the various columns (or cells) in the .csv file to variables in R by their column headings
p0 <- d$p0
m0 <- d$m0
pstar <- d$pstar
mstar <- d$mstar
mortparms <- d$mortparms
fertparm <- d$fertparm[1]
rho <- d$rho
phi <- d$phi
mortscale <- d$mortscale[1]
fertscale <- d$fertscale[1]
N.AC <- length(p0)
Y <- d$Y[1]*0.9
Y.CV <- d$Y.CV[1]
Yac <- 0.3
Am <- d$Am[1]
J <- d$J[1]
h <- d$h[1]
N.guess <- d$N.guess[1]
span <- d$span[1]

#################################################################################################### Set Pop2 parms
# Assign the various columns (or cells) in the .csv file to variables in R by their column headings
p0b <- d$p0
m0b <- d$m0
pstarb <- d$pstar
mstarb <- d$mstar
mortparmsb <- d$mortparms
fertparmb <- d$fertparm[1]
rhob <- d$rho
phib <- d$phi
mortscaleb <- d$mortscale[1]
fertscaleb <- d$fertscale[1]
N.ACb <- length(p0)
Yb <- d$Y[1]
Y.CVb <- d$Y.CV[1]
Yacb <- 0.3
Amb <- d$Am[1]
Jb <- d$J[1]
hb <- d$h[1]
Nb.guess <- d$N.guess[1]

############# IDD Parms
P1 <- P2 <- Q <- Pi1 <- Pi2 <- M1 <- M2 <- c()

c <- 1
g <- 0.1
y <- 1

m <- 1
v <- 2

R <- 100

runlength<-1

###################################################################################################################################
####################################################################################################### Load Functions and Packages
    
# load.sources.r must be edited to point to the right directory
source(paste0(SS,"/load.sources.r"))
load.sources("Sources.r")

require(pracma)

###################################################################################################################################
################################################################################################# Create list to store results

# this is a list
RUN <- vector("list",5)

# Equivalent to Nstoch, nstoch and Estoch from Matlab code, or (1) total pop size, (2) pop vector and (3) E, each by time
K <- list(Ystoch=matrix(NA,ncol=span,nrow=runlength),NstochNM=matrix(NA,ncol=span,nrow=runlength),nstochNM=array(NA,dim=c(N.AC,span,runlength)), Nstoch=matrix(NA,ncol=span,nrow=runlength),nstoch=array(NA,dim=c(N.AC,span,runlength)),Estoch=matrix(NA,ncol=span,nrow=runlength),M=array(NA,dim=c(N.AC,span,runlength)))

Kb <- list(Ystoch=matrix(NA,ncol=span,nrow=runlength),NstochNM=matrix(NA,ncol=span,nrow=runlength),nstochNM=array(NA,dim=c(N.AC,span,runlength)),Nstoch=matrix(NA,ncol=span,nrow=runlength),nstoch=array(NA,dim=c(N.AC,span,runlength)),Estoch=matrix(NA,ncol=span,nrow=runlength),M=array(NA,dim=c(N.AC,span,runlength)))


###################################################################################################################################
####################################################################################################### Find Equilibrium Values

# Assign the output of the F.Ehat function to the multidimensonal array "G"
G <- F.Ehat(p0,m0,pstar,mstar,mortparms,fertparm,rho,phi,fertscale,mortscale,N.AC)

# find Nhat, given Ehat and population structure at equilibrium
scrap.Nhat <- fzero(function(Z) ((Y*Am*(1-exp(-h*G$phi1*Z)))/J/G$rho1/Z-G$Ehat), N.guess)

# define Nhat or give error message if the population is negative
Nhat <- ifelse(scrap.Nhat$x<0,(stop("Not a viable equilibrium population")),scrap.Nhat$x)

# nhat is the vector of equilibrium population by age
nhat <- Nhat*G$uhat;

e0hat <- demogcalc(G$phat,G$mhat)[[1]]
e10hat <- demogcalc(G$phat,G$mhat)[[2]]
TFRhat <- demogcalc(G$phat,G$mhat)[[3]]
NRRhat <- demogcalc(G$phat,G$mhat)[[4]]
muChat <- demogcalc(G$phat,G$mhat)[[5]]


###################################################################################################################################
###################################################################################################################################
####################################################################################################### Begin Time Series Analysis

# do multiple runs of each yield modification
for (r in 1:1){
   
  # loop through each yield scenario
  for (j in 1:runlength){
    tryCatch({
    ############################################################################
    ################################################# Set Inital Conditions
    
    #################### Population 1
      # use a seed to allow replication of stochastic series  
        set.seed(666+(r*500))
      # Determine "baseline" (full-food) projection matrix
      A.baseline <- quickleslie(m0,p0,N.AC)
      
      # Define initial population size and vector
      K$Nstoch[j,1] <- Nhat/50
      K$nstoch[,1,j] <- G$uhat*K$Nstoch[j,1]
      
      # Draw a yield from gamma distrib w mean=Y and CV=Y.CV for 1st entry only here
      K$Ystoch[j,1] <- rgamma(1,shape=(1/(Y.CV^2)),scale=(Y*(Y.CV^2)));
      
      #################### Population 2
      # use a seed to allow replication of stochastic series  
        set.seed(666+(r*500))
      # Determine "baseline" (full-food) projection matrix
      A.baselineb <- quickleslie(m0b,p0b,N.ACb)
      
      # Define initial population size and vector
      Kb$Nstoch[j,1] <- Nhat/50
      Kb$nstoch[,1,j] <- G$uhat*Kb$Nstoch[j,1]
      
      # Draw a yield from gamma distrib w mean=Y and CV=Y.CV for 1st entry only here
      Kb$Ystoch[j,1] <- rgamma(1,shape=(1/(Y.CVb^2)),scale=(Yb*(Y.CVb^2)));
      
      #################### IDD

      
    ######################################################################################################################################
    ################################################################################################################### Iterate over years
    
      for(i in 2:span){
      ###################################################################################################### Model Peaseants of Population 1
       set.seed(666 + i + r*500) 
        # Draw a yield from gamma distrib w mean=Y and CV=Y.CV
        Ybar <- Y*(1-Yac)
        Ystar <- Ybar + Yac*K$Ystoch[j,(i-1)]
        K$Ystoch[j,i] <- rgamma(1,shape=(1/(Y.CV^2)),scale=(Ystar*(Y.CV^2)));
        
        # Calculate fraction of land currently in cultivation (F.scrap here), then E(t)
        F.scrap <- (1-exp(-h*sum(K$nstoch[,(i-1),j]*phi)))
        K$Estoch[j,(i-1)] <- (K$Ystoch[j,i]*Am*F.scrap)/(J*sum(K$nstoch[,(i-1),j]*rho))
        
        # Given Estoch, calculate the projection matrix and step population forward
        if( is.na(K$Estoch[j,(i-1)])){
          
          stop(paste(c("The 1st population in run ", j, " has crashed to extinction and is not viable.")))
        }
        
        if(K$Estoch[j,(i-1)]>=1){
          
          # If Estoch >= 1, just use baseline projection matrix
          K$nstoch[,i,j] <- A.baseline%*%K$nstoch[,(i-1),j]
          K$Nstoch[j,i] <- sum(K$nstoch[,i,j])
        }
        else{
          # If Estoch < 1, use Gamma cdf to calc food-limited vital rates for projection matrix
          # calc survival and fertility at given value of E (expressed as Z here)
          scrap.p2 <- pstar*pgamma(K$Estoch[j,(i-1)], shape=mortparms, scale=mortscale)
          scrap.m2 <- mstar*pgamma(K$Estoch[j,(i-1)], shape=fertparm, scale=fertscale)
          A.stoch <- quickleslie(scrap.m2,scrap.p2,N.AC)
          
          K$nstoch[,i,j] <- A.stoch%*%K$nstoch[,(i-1),j]
           K$Nstoch[j,i] <- sum(K$nstoch[,i,j])
        } # End if else
        
        ###################################################################################################### Model Peaseants of Population 2
        set.seed(666 + i + r*500) 
        # Draw a yield from gamma distrib w mean=Y and CV=Y.CV
        Ybarb <- Yb*(1-Yacb)
        Ystarb <- Ybarb + Yacb*Kb$Ystoch[j,(i-1)]
        Kb$Ystoch[j,i] <- rgamma(1,shape=(1/(Y.CVb^2)),scale=(Ystarb*(Y.CVb^2)));
        
        # Calculate fraction of land currently in cultivation (F.scrap here), then E(t)
        F.scrapb <- (1-exp(-hb*sum(Kb$nstoch[,(i-1),j]*phib)))
        Kb$Estoch[j,(i-1)] <- (Kb$Ystoch[j,i]*Amb*F.scrapb)/(Jb*sum(Kb$nstoch[,(i-1),j]*rhob))
        
        # Given Estoch, calculate the projection matrix and step population forward
        if( is.na(Kb$Estoch[j,(i-1)])){
          
          stop(paste(c("The 2nd population in run ", j, " has crashed to extinction and is not viable.")))
        }
        
        if(Kb$Estoch[j,(i-1)]>=1){
          
          # If Estoch >= 1, just use baseline projection matrix
          Kb$nstoch[,i,j] <- A.baselineb%*%Kb$nstoch[,(i-1),j]
          Kb$Nstoch[j,i] <- sum(Kb$nstoch[,i,j])
        }
        else{
          # If Estoch < 1, use Gamma cdf to calc food-limited vital rates for projection matrix
          # calc survival and fertility at given value of E (expressed as Z here)
          scrap.p2b <- pstarb*pgamma(Kb$Estoch[j,(i-1)], shape=mortparmsb, scale=mortscaleb)
          scrap.m2b <- mstarb*pgamma(Kb$Estoch[j,(i-1)], shape=fertparmb, scale=fertscaleb)
          A.stochb <- quickleslie(scrap.m2b,scrap.p2b,N.ACb)
          
          Kb$nstoch[,i,j] <- A.stochb%*%Kb$nstoch[,(i-1),j]
          Kb$Nstoch[j,i] <- sum(Kb$nstoch[,i,j])
        } # End if else

      ###############################################################################################################################
      ################################################################################# Model IDD Dynamics        
         Pi1[i] <- max(0, K$Estoch[j,i-1] )
         Pi2[i] <- max(0, Kb$Estoch[j,i-1] )
         
         P1 <- K$nstochNM[,i,j] <- K$nstoch[,i,j]
         P2 <- Kb$nstochNM[,i,j] <- Kb$nstoch[,i,j]
         
         K$NstochNM[j,i] <- sum(K$nstoch[,i,j])
         Kb$NstochNM[j,i] <- sum(Kb$nstoch[,i,j])
         
         #### Allow Migration
         for(n in 1:N.AC){
         P1[n] <- P1[n]  -  (m*P1[n]*((Pi2[i]^v)/((c*Pi1[i])^v+((Pi2[i])^v))))  +  (m*P2[n]*(((c*Pi1[i])^v)/((c*Pi1[i])^v+((Pi2[i])^v)))) 

         P2[n] <- P2[n]  +  (m*P1[n]*((Pi2[i]^v)/((c*Pi1[i])^v+((Pi2[i])^v))))  -  (m*P2[n]*(((c*Pi1[i])^v)/((c*Pi1[i])^v+((Pi2[i])^v))))

         M1[n]  <-           (m*P1[n]*((Pi2[i]^v)/((c*Pi1[i])^v+((Pi2[i])^v))))  -  (m*P2[n]*(((c*Pi1[i])^v)/((c*Pi1[i])^v+((Pi2[i])^v))))
         M2[n]  <-           (m*P2[n]*(((c*Pi1[i])^v)/((c*Pi1[i])^v+((Pi2[i])^v)))) - (m*P1[n]*((Pi2[i]^v)/((c*Pi1[i])^v+((Pi2[i])^v))))
         
         if( P1[n]<M1[n] | P1[n]<M1[n]){
          stop(paste(c("The migration rate parameter is too high ", "more people left the population than live in the population.")))
        }
         }
         
         K$nstoch[,i,j] <- P1
         Kb$nstoch[,i,j] <- P2
         
         K$M[,i,j] <- M1
         Kb$M[,i,j] <- M2
         
         K$Nstoch[j,i] <- sum(K$nstoch[,i,j])
         Kb$Nstoch[j,i] <- sum(Kb$nstoch[,i,j])
  
      ###############################################################################################################################
      ###################################################################################################################  End Model    
      } # End span loop
    } # End error catching function
    , error=function(e){cat("Warning:",conditionMessage(e), "\n")})
  } # End runlength loop
  RUN[[r]] <- K
} # End RUN (seed) loop


###################################################################################################################################
####################################################################################################### Plotting and Tables
  plot(Kb$Estoch[j,]-K$Estoch[j,])
# Plot time series of population size under a variety of food yield ratios.
# ONLY plots last loop through r. Use data stored in RUN to compare the various Ks
plot(K$Nstoch[1,], type="l", ylim=c(0, max(K$Nstoch[,],na.rm=T)),lwd=1,col="red")
lines(Kb$Nstoch[1,], col="black")

windows()
mig<-c()
for(i in 1:length(K$Nstoch))
mig[i]<-sum(K$M[,i,1])/K$Nstoch[i]

plot(mig,type="l")





