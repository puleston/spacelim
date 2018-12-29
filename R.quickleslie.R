
quickleslie <- function(mvec,pvec,N.AC){

# initialize the Leslie matrix w zeros
scrap.matrix <- matrix(0,nrow=N.AC,ncol=N.AC)

# fill off-diagonal w survival probabilities
for(i in 1:(N.AC-1)){
	scrap.matrix[i+1,i]<-pvec[i]
	}

# put fert probabilities in first row of matrix
for(i in 1:N.AC){
	scrap.matrix[1,i]<-mvec[i]
	}
	
scrap.matrix
 
 	}