
################################################################################## R model
# New function F.Ehat

# p0 - vector - survival prob from age class to next at full food (only use when E>=1)
# m0 - vector - fertility prob from age class to next at full food (only use when E>=1)
# pstar - vector - baseline survival prob for use w E < 1, they make gamma cdf give correct answer when E=1 
# mstar - vector - baseline fertility prob for use w E < 1, they make gamma cdf give correct answer when E=1 
# mortparms - vector - alpha of gamma cdf, responsivness of survival on food
# fertparm - scalar - alpha of gamma cdf, responsiveness of fertility to food
# rho - vector - consumption for each age class, relative to highest consuming age class
# phi - vector - worker effectivness, relative to most effective age class
# N.AC - scalar, number of age classes

F.Ehat <-function(p0,m0,pstar,mstar,mortparms,fertparm,rho,phi,fertscale,mortscale, N.AC){

#Run Calculations on Input
 require(pracma)
 
# Run function eqmmatrixfun_cp, as defined before
 scrap.Ehat <- fzero(function(Z) (eqmmatrixfun_cp(Z,mstar,pstar,fertparm,mortparms,fertscale,mortscale,N.AC)-1), c(0,1))
 
 Ehat<-scrap.Ehat$x


# fhat is equilibrium survival multipliers as a function of food availability
 fhats <-pgamma(Ehat, shape=mortparms, scale=mortscale)
# fhats[N.AC]<-0

# ghat is equilibrium fertility multipliers as a function of food availability
 ghat <-pgamma(Ehat, shape=fertparm, scale=fertscale)

# qhat is 1 minus elementwise product; probability of *not* surviving from one age class to the next at equilibrium;
 qhat <- (1-(pstar*fhats))

# phat is vector of equilibrium survival probabilities; p(Ehat)
 phat <- pstar*fhats
 
# mhat is vector of equilibrium fertility probabilities; m(Ehat)
 mhat<- mstar*ghat

# lhat is cumulative survival probs to age x evaluated at equilibrium
 lhat <-c() #initialize the variable
 scrap1 <- cumprod(phat); #calculate the cumulative products of the surival probability vector
 lhat[1] <- 1 #1st entry of the lhat vector is 1, indicating that probability of surviving to the first age class is assumed to be 1, given that there is a live birth, second entry is p(1), third is p(1)*p(2), etc.
 lhat[2:N.AC]<- scrap1[1:(N.AC-1)] #assemble the vector

# Ahat is a leslie matrix using other custom function
 Ahat <- quickleslie(mhat,phat,N.AC); 

# Calc all eigen vals and vectors
 scrap2 <- eigen(Ahat)

 # Calc dom eigenvalue by finding the one of the greatest magnitude (absolute value)
 scrap.loc.eig <- which(abs(scrap2$values)==max(abs(scrap2$values)))
 scrap.loc.eig <- scrap.loc.eig[1]
 
 # Calc dom eigenvector by grabbing the vector that corresponds to the index of the dominant eigenvalue
 dom.vec1 <- scrap2$vec[,scrap.loc.eig]
 # equilibrium age structure is proportional to the dominant eigenvector, called uhat
 uhat <- abs(dom.vec1/(sum(dom.vec1)))

# phi1 is dot product of row vectors with uhat, represents the equilibrium per capita relative work effectiveness; rho1 is the same, but for consumption; both are scalars
 phi1 <- phi%*%uhat
 rho1 <- rho%*%uhat
 
# Return Output
list(Ehat=Ehat,qhat=qhat,phat=phat,mhat=mhat,lhat=lhat,uhat=uhat,Ahat=Ahat,phi1=phi1,rho1=rho1)
}



