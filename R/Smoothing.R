
spatialSmoothing <- function(spatSmooth, outbreakData2, lambda) {
	if (spatSmooth)
	{
		dim = dim(outbreakData2$data)   # dimension of the space (Y x X x T)
		T = dim[3]
		
		Obs3 <- array(rep(0,prod(dim)), dim)
		Obs3.mu <- array(rep(0,prod(dim)), dim)
		
		for (t in 1:T)
		{
			# data matrix containing observed values:
			data <- outbreakData2$data[,,t]
			data <- as.matrix(data)             # the data matrix to be partitioned
			Obs3[,,t] <- spatial.smooth(data,lambda)
			
			# corresponding matrix containing expected values:
			data.mu <- outbreakData2$means[,,t]
			data.mu <- as.matrix(data.mu)
			Obs3.mu[,,t] <- spatial.smooth(data.mu,lambda)
		}
	} else {
		Obs3 <- outbreakData2$data
		Obs3.mu <- outbreakData2$means
	}
	
	return (list(data=Obs3, means=Obs3.mu))
}

spatial.smooth <-function (Yt, lambda){
	nrow = dim(Yt)[1]
	ncol = dim(Yt)[2]
	
	SA <- matrix(0,nrow,nrow)
	for(i in 1:nrow) {
		for(j in 1:nrow) {
			SA[i,j] <- lambda * (1-lambda)^{abs(i-j)}
		}
	}
	v1A <- rep(1, nrow)
	
	SB <- matrix(0,ncol,ncol)
	for(i in 1:ncol) {
		for(j in 1:ncol) {
			SB[i,j] <- lambda * (1-lambda)^{abs(i-j)}
		}
	}
	v1B <- rep(1, ncol)
	
	totA <- SA %*% v1A
	DA <- diag(c(totA))
	iDA <- diag(1/c(totA))  # DA^(-1)
	
	totB <- SB %*% v1B
	DB <- diag(c(totB))
	iDB <- diag(1/c(totB))  # DB^(-1)
	
	Yt2 <- iDA %*% SA %*% Yt %*% SB %*% iDB
	
	return (Yt2)
}

temporalSmoothing <- function(tempSmooth, outbreakData, alpha) {
	if (tempSmooth)	{
		data <- temporal.smooth2(outbreakData$data, alpha, outbreakData$data0)
		
		means <- temporal.smooth2(outbreakData$means, alpha, outbreakData$means0)
	} else {
		data <- outbreakData$data
		means <- outbreakData$means
	}

	return (list(data=data, means=means))
}

temporal.smooth2 <-function (Yt, alpha, means0) {
	# dimension of the space (Y x X x T)
	Y = dim(Yt)[1]
	X = dim(Yt)[2]
	T = dim(Yt)[3]
	
	dim = c(Y,X,T)
	Zt <- array(rep(0,prod(dim)), dim)
	
	for (i in 1:Y) {
		for (j in 1:X) {
			prevZ = means0[i,j]
			for (t in 1:T) {
				Zt[i,j,t] = alpha * Yt[i,j,t] + (1-alpha) * prevZ;
				prevZ = Zt[i,j,t]
			}
		}
	}
	return (Zt)
}

