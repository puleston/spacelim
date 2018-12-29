
eqmmatrixfun_cp <- function(Z,mstar,pstar,fertparm,mortparms,fertscale,mortscale,N.AC){

# calc survival and fertility at given value of E (expressed as Z here)	
scrap.p <- pstar*pgamma(Z, shape=mortparms, scale=mortscale)
scrap.m <- mstar*pgamma(Z, shape=fertparm, scale=fertscale)

# initialize the Leslie matrix w zeros
scrap.matrix <- matrix(0,nrow=N.AC,ncol=N.AC)

# fill off-diagonal w survival probabilities
for(i in 1:(N.AC-1)){
	scrap.matrix[i+1,i]<-scrap.p[i]
	}

# put fert probabilities in first row of matrix
for(i in 1:N.AC){
	scrap.matrix[1,i]<-scrap.m[i]
	}
	
# Calc all eigen vals and vectors
 scrap.eig <- eigen(scrap.matrix)

 # Calc dom eigenvalue by finding the one of the greatest magnitude (absolute value)
 scrap.lambda <- scrap.eig$values[which(abs(scrap.eig$values)==max(abs(scrap.eig$values)))]
 # Take the first value if there are two of identical magnitude
 scrap.lambda <- scrap.lambda[1]
 
 abs(scrap.lambda)
 
 	}