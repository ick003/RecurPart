# This file works with Outbreak_Detection4.R, Outbreak_Detection5.R & Plot_AEGISS.R

### Generate time dependance in Poisson Counts
# inputs:
#   mu = average ln Poisson mean
#   bweekday = proportions of influence for each day of week, sum to 0
#   bcos, bsin: radius of influence
#   nob = number of days to be generated
#   blags = extra vector of diminishing influence, e.g. for 4 days
#   k = min mean value allowed in Poisson distribution
#
# returning ncount = a vector of poisson counts for nob days
#
transitionalPoisson<-function(mu,bweekday,bcos,bsin,blags,nob,k){
  nlag<-length(blags)   # nlag = number of diminishing days that can influence the Poisson mean
  lcounts<-rep(0,nlag)  # lcounts: e.g. latest 4 days of Poisson counts influencing Poisson mean
  
  dw<-1   # dw = day of week (1 to 7)
  for (i in 1:nob) {
  	# temp = ln of mean value of the next Poisson distribution
    temp <- mu+bweekday[dw]+bcos*cos(2*pi*i/365.25)+bsin*sin(2*pi*i/365.25)+
    	sum(blags*log(lcounts+1))
    
    temp[temp<k]<-k
    
    # ncount = current generated Poisson count
    if(i==1)
    	ncount<-rpois(1,exp(temp))
    else
    	ncount<-c(ncount,rpois(1,exp(temp)))
    
    lcounts[2:nlag]<-lcounts[1:(nlag-1)] # shift lcounts right by 1 number
    lcounts[1]<-ncount[i]     # store latest generated Poisson count into 1st value of lcounts
    
    dw<-dw+1  # cycle day of week
    if(dw>7)dw<-1
  }
  
  ncount
}

#### Applying the smoothing matrix
# inputs:
#   x = temporally smoothed matrix Yt
#   F40 = DA^(-1) x SA
#
# returning temp = DA^(-1) x SA x Yt- x SB x DA^(-1) = Yt~
#   i.e. spatially smoothed EWMA counts
spatial.ewma.smooth<-function(x, F40) {
  temp<-F40%*%x
  temp<-temp%*%t(F40)
  
  return(temp)
}

temporal.smooth <-function (Yt,alpha,Z0) {
	X = dim(Yt)[1]
	Y = dim(Yt)[2]
	T = dim(Yt)[3]
	
	dim = c(X,Y,T)
	Zt <- array(rep(0,prod(dim)), dim)
	
	for (i in 1:X) {
		for (j in 1:Y) {
			prevZ = Z0
			for (t in 1:T) {
				Zt[i,j,t] = alpha * Yt[i,j,t] + (1-alpha) * prevZ;
				prevZ = Zt[i,j,t]
			}
		}
	}
	return (Zt)
}

temporal.smooth2 <-function (Yt,alpha,means) {
	X = dim(Yt)[1]
	Y = dim(Yt)[2]
	T = dim(Yt)[3]
	
	dim = c(X,Y,T)
	Zt <- array(rep(0,prod(dim)), dim)
	
	for (i in 1:X) {
		for (j in 1:Y) {
			prevZ = means[i,j]
			for (t in 1:T) {
				Zt[i,j,t] = alpha * Yt[i,j,t] + (1-alpha) * prevZ;
				prevZ = Zt[i,j,t]
			}
		}
	}
	return (Zt)
}

spatial.smooth <-function (Yt,lambda){
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

#### Generate Time-dependant Poisson Counts matrices
# inputs:
#   N = dimension for x- and y-coordinates => a square
#   T = number of days to generate Poisson counts for each x and y
#
# returning counts = Poisson counts in 3D (t:1..T, x:1..N, y:1..N)
#
gen.Poisson <- function(N,T){
	counts<-array(0,c(T,N,N))
	
	for(i in 1:N) {
	  for(j in 1:N) {
	    counts[,i,j] <- transitionalPoisson(mu=-0.123,bweekday=c(rep(-0.1,5),0.2,0.3)*0.25,
	    							    bcos=0.1,bsin=-0.1,blags=c(0.1,0.06,0.02,0.01),nob=T,k=-4)
	  }
	}
	
	return(counts)
}

#### Generate Spatial dependance in Poisson Counts
# inputs:
#   counts = Poisson counts in 3D (t:1..T, x:1..N, y:1..N)
#   T = number of days to generate Poisson counts for each x and y
#   N = dimension for x- and y-coordinates => a square
#
# returning sp.counts = Poisson counts in 3D (t:1..T, x:1..N, y:1..N) with spatial dependence
#
gen.Poisson.sp <- function(counts,T,N){
	# sp.counts = Poisson counts in 3D (t:1..T, x:1..N, y:1..N) with spatial dependance / influence
	sp.counts<-array(0,c(T,N,N))
	
	for(i in 1:(N-1)) {
	  for(j in 1:(N-1)) {
	    if (i==1) {
	    	if(j==1)
	    		sp.counts[,i,j]<-floor(counts[,i,j]+
	    			0.3*(counts[,i+1,j]+counts[,i,j+1])+
	    			0.1*counts[,i+1,j+1]) 
	    	else 
	    		sp.counts[,i,j]<-floor(counts[,i,j]+
	    			0.3*(counts[,i+1,j]+counts[,i,j+1]+counts[,i,j-1])+
	    			0.1*(counts[,i+1,j+1]+counts[,i,j-1]))  #? 0.1*(counts[,i+1,j+1]+counts[,i+1,j-1]))
	    } else {
	      if (j==1)
	      	sp.counts[,i,j]<-floor(counts[,i,j]+
	      		0.3*(counts[,i+1,j]+counts[,i,j+1]+counts[,i-1,j])+
	      		0.1*(counts[,i+1,j+1]+counts[,i-1,j+1]))
	     	else 
	        sp.counts[,i,j]<-floor(counts[,i,j]+
	        	0.3*(counts[,i+1,j]+counts[,i,j+1]+counts[,i,j-1]+counts[,i-1,j])+
	        	0.1*(counts[,i-1,j+1]+counts[,i+1,j+1]+counts[,i+1,j-1]+counts[,i-1,j-1]))
	    } 
	  }
	}
	
	i<-N
	for(j in 1:(N-1)){
	  if(j==1) 
	  	sp.counts[,i,j] <- floor(counts[,i,j]+
	  													 	0.3*(counts[,i-1,j]+counts[,i,j+1])+
	  													 	0.1*counts[,i-1,j+1])
	  else 
	  	sp.counts[,i,j] <- floor(counts[,i,j]+
	  													 	0.3*(counts[,i-1,j]+counts[,i,j+1]+counts[,i,j-1])+
	  													 	0.1*(counts[,i-1,j+1]+counts[,i-1,j-1]))
	}
	
	j<-N
	for(i in 1:(N-1)){
	  if(i==1)
	  	sp.counts[,i,j] <- floor(counts[,i,j]+
	  													 	0.3*(counts[,i+1,j]+counts[,i,j-1])+
	  													 	0.1*counts[,i+1,j-1]) 
	  else 
	  	sp.counts[,i,j] <- floor(counts[,i,j]+
	  													 	0.3*(counts[,i-1,j]+counts[,i,j-1]+counts[,i+1,j])+
	  													 	0.1*(counts[,i+1,j-1]+counts[,i-1,j-1]))
	}
	
	i<-N;j<-N
	sp.counts[,i,j] <- floor(counts[,i,j]+
													 	0.3*(counts[,i-1,j]+counts[,i,j-1])+
													 	0.1*counts[,i-1,j-1])
	
	return(sp.counts)
}

#### Generate Outbreak in Poisson Counts
# inputs:
#   sp.counts = Poisson counts in 3D (t:1..T, x:1..N, y:1..N) with spatial dependence
#   xyt.outbreak = outbreak size in 3D
#
# returning sp.counts = modified Poisson counts in 3D with generated outbreaks in volume
#   defined by xyt.outbreak
#
gen.Poisson.outbreak <- function(sp.counts,xyt.outbreak){
  xoutbreak = xyt.outbreak$x
  youtbreak = xyt.outbreak$y
  t.outbreak = xyt.outbreak$t
  
  sp.counts[t.outbreak,xoutbreak,youtbreak] <- 
  	sp.counts[t.outbreak,xoutbreak,youtbreak]+
  	  array(rpois(length(t.outbreak)*length(xoutbreak)*length(youtbreak),delta),
  	  c(length(t.outbreak),length(xoutbreak),length(youtbreak)))
  
	return(sp.counts)
}

#### Generate EWMA Temporal Smoothing 
# inputs:
#   sp.counts = Poisson counts in 3D with spatial dependence and outbreaks
#   weight = exponential weights for days 1..T
#   time = days in 1..T
#   dw = days of weeks (1..7) in 1:T
#   RL = run length for delayed processing
#   muo = mean of Poisson counts
#   T0 = offset days to start for Poisson fitting
#
# returning ewma in 3D (1, 1..N, 1..N) for observed counts with temporal smoothing
#   ewma.mu for expected counts with temporal smoothing
temp.smooth.ewma <- function(sp.counts, weight, time, dw, RL, muo, T0) {
	N <- length(sp.counts[1,1,])
	
	ewma<-array(0,c(T-T0,N,N))
	ewma.mu<-array(0,c(T-T0,N,N))
	
	for(i in 1:N) {
	  for(j in 1:N) {
	    times<-time[c(5:T0)]   # times = 5..T0
	    wds<-dw[c(5:T0)]       # wds = days of week 5..T0
	    
	    lag4<-log(sp.counts[1:(T0-4),i,j]+1)  # lag 4 days from T0
	    lag3<-log(sp.counts[2:(T0-3),i,j]+1)  # lag 3 days from T0
	    lag2<-log(sp.counts[3:(T0-2),i,j]+1)  # lag 2 days from T0
	    lag1<-log(sp.counts[4:(T0-1),i,j]+1)  # lag 1 day from T0
	    
	    # fitting a Poisson distribution to the Poisson counts RL+4..T0-1+RL
	    glm.4 <- try(glm(sp.counts[(RL+4):(T0-1+RL),i,j]~1 + factor(wds) + 
	    								 	lag1 + lag2 + lag3 + lag4 + 
	    								 	cos(times*2*pi/365.25)+sin(times*2*pi/365.25), 
	    								  weights=weight[1:(T0-4)], family=poisson()))
	    
	    ewma[1,i,j] <- 0.1*sp.counts[1,i,j]+0.9*muo
	    ewma.mu[1,i,j] <- 0.1*predict(glm.4, type="response")[T0-4]+0.9*muo
	  }
	}
	
	return(list(ewma = ewma, ewma.mu = ewma.mu))
}

# calculate sums in (row) X-, (column) Y- and (time) T-dimensions respectively:
# (not cummulative sums)
dimSum <- function(scounts) {
	nX <- dim(scounts)[1]
	nY <- dim(scounts)[2]
	nT <- dim(scounts)[3]

	sumX <- rep(0,nX)
	index <- 0
	for (i in 1:nX)
	{
		index <- index+1
		for (j in 1:nY)
		{
			for (k in 1:nT)
			{
				sumX[index] = sumX[index] + scounts[i,j,k]
			}
		}
	}

	sumY <- rep(0,nY)
	index <- 0
	for (j in 1:nY)
	{
		index <- index+1
		for (i in 1:nX)
		{
			for (k in 1:nT)
			{
				sumY[index] = sumY[index] + scounts[i,j,k]
			}
		}
	}
	
	sumT <- rep(0,nT)
	index <- 0
	for (k in 1:nT)
	{
		index <- index+1
		for (i in 1:nX)
		{
			for (j in 1:nY)
			{
				sumT[index] = sumT[index] + scounts[i,j,k]
			}
		}
	}
	
	msum <- list(X=sumX,Y=sumY,T=sumT)
	return (msum)
}
	
# inputs:
#       scounts: the current 3D array O2 to be analysed for row/col/tim partitioning
#       smus: corresponding 3D array with expected counts (same dimension as scounts)
#       gen.dim: $nX: initial row and last row indices
#                $nY: initial col and last col indices
#                $nT: initial tim and last tim indices
#
# choosing the best row/col/tim partition A for 1 iteration, and returning part: 
#       type.part = 1 if choosing row partition
#                   2 if choosing col partition
#                   3 if choosing tim partition
#       partition = [1] = initial index for row/col/tim partition A
#                   [2] = last index for row/col/tim partition
#       gen.dim: $nX / $nY / $nT updated with the chosen partition A
#
# partitioning by z-scores only (i.e. not by p-values):
recursive.partition3D <- function(scounts,smus,gen.dim) {
	# yr, yrmu: row sums
	# yc, ycmu: col sums
	# yt, ytmu: tim sums
	msum <- dimSum(scounts)
	yr <- msum$X; yc <- msum$Y; yt <- msum$T
	
	msum.mu <- dimSum(smus)
	yrmu <- msum.mu$X; ycmu <- msum.mu$Y; ytmu <- msum.mu$T
	
	nX <- length(yr)
	nY <- length(yc)
	nT <- length(yt)
	
	# pr.up: z-scores in upward row direction, pr.down: z-scores in downward row direction
	# pc.up: z-scores in rightward col direction, pc.down: z-scores in leftward col direction
	pr.up   <- 2*(sqrt(cumsum(yr))-sqrt(cumsum(yrmu)))
	pr.down <- 2*(sqrt(cumsum(yr[nX:1]))-sqrt(cumsum(yrmu[nX:1])))
	
	pc.up   <- 2*(sqrt(cumsum(yc))-sqrt(cumsum(ycmu)))
	pc.down <- 2*(sqrt(cumsum(yc[nY:1]))-sqrt(cumsum(ycmu[nY:1])))
	
	pt.up   <- 2*(sqrt(cumsum(yt))-sqrt(cumsum(ytmu)))
	pt.down <- 2*(sqrt(cumsum(yt[nT:1]))-sqrt(cumsum(ytmu[nT:1])))
	
	# make last z-score = min. z-score - 1 in order to eliminate the case of no partitioning:
	zmin = min(c(pr.up, pr.down, pc.up, pc.down, pt.up, pt.down))
	
	pr.up[nX] = zmin-1; pr.down[nX] = zmin-1
	pc.up[nY] = zmin-1; pc.down[nY] = zmin-1
	pt.up[nT] = zmin-1; pt.down[nT] = zmin-1
	
	# tem.X: max z-scores for row partition
	# tem.Y: max z-scores for col partition
	# tem.T: max z-scores for tim partition
	tem.X <- c(max(pr.up),max(pr.down))
	tem.Y <- c(max(pc.up),max(pc.down))
	tem.T <- c(max(pt.up),max(pt.down))
	
	# Look for the max z-score for row sums - summing over rows up and then down
	# p1: index for max z-score for row partition in upward direction
	# p2: index for max z-score for row partition in downward direction
	p1 <- c(1:nX)[pr.up==tem.X[1]]
	p2 <- c(nX:1)[pr.down==tem.X[2]]
	
	# Look for the max z-score for col sums - summing over cols up and then down
	# p3: index for max z-score for col partition in upward direction
	# p4: index for max z-score for col partition in downward direction
	p3<-c(1:nY)[pc.up==tem.Y[1]]
	p4<-c(nY:1)[pc.down==tem.Y[2]]
	
	# Look for the max z-score for tim sums - summing over times up and then down
	# p5: index for max z-score for tim partition in upward direction
	# p6: index for max z-score for tim partition in downward direction
	p5<-c(1:nT)[pt.up==tem.T[1]]
	p6<-c(nT:1)[pt.down==tem.T[2]]
	
	# sometimes, there are multiple max z-scores
	if (length(p1)>1) p1 <- p1[1]
	if (length(p2)>1) p2 <- p2[1]
	if (length(p3)>1) p3 <- p3[1]
	if (length(p4)>1) p4 <- p4[1]
	if (length(p5)>1) p5 <- p5[1]
	if (length(p6)>1) p6 <- p6[1]
	
	# Recording all potential "best" partitions
	# potential2: 1st row index, row index for max upward z-score, 
	#                            row index for max downward z-score, last row index,
	#             1st col index, col index for max leftward z-score,
	#                            col index for max rightward z-score, last col index.
	potential2 <- matrix(c(gen.dim$nX[1],gen.dim$nX[1]+p1-1,
												 gen.dim$nX[1]+p2-1,gen.dim$nX[2],
												 gen.dim$nY[1],gen.dim$nY[1]+p3-1,
												 gen.dim$nY[1]+p4-1,gen.dim$nY[2],
											   gen.dim$nT[1],gen.dim$nT[1]+p5-1,
											   gen.dim$nT[1]+p6-1,gen.dim$nT[2]),
												 6,2,byrow=T)
	
	# Check whether no partitioning is necessary of the rows
	trmax <- max(tem.X)  # overall max z-score in row partition
	row1 <- 1 # row1 = 1 to signify best row partition in upward direction
	if (gen.dim$nX[1]!=gen.dim$nX[2]) {
		if (tem.X[1]!=trmax) {
			row1 <- 2 # row1 = 2 to signify best row partition in downward direction
		}
	} 
	
	# Check whether no partitioning is necessary of the columns
	tcmax <- max(tem.Y)  # overall max z-score in col partition
	row2 <- 3 # row2 = 3 to signify best col partition in upward direction
	if (gen.dim$nY[1]!=gen.dim$nY[2]) {
		if (tem.Y[1]!=tcmax) {
			row2 <- 4 # row2 = 4 to signify best col partition in downward direction
		}
	} 

	# Check whether no partitioning is necessary of the times
	ttmax <- max(tem.T)  # overall max z-score in tim partition
	row3 <- 5 # row3 = 5 to signify best tim partition in upward direction
	if (gen.dim$nT[1]!=gen.dim$nT[2]) {
		if (tem.T[1]!=ttmax) {
			row3 <- 6 # row3 = 6 to signify best tim partition in downward direction
		}
	}
	
	# tmax: overall max normalized z-score for both row & col partitions
	# ind: 1 to signify best row partition in upward direction
	#      2 to signify best row partition in downward direction
	#      3 to signify best col partition in leftward direction
	#      4 to signify best col partition in rightward direction
	#      5 to signify best tim partition in upward direction
	#      6 to signify best tim partition in downward direction
	# col.part: 1 if choosing row partition
	#           2 if choosing col partition
	#           3 if choosing tim partition
	tmax = max(tem.X, tem.Y, tem.T)
	if(trmax==tmax) {
		ind<-row1
		col.part<-1
		if (ind==1) {
			A_value = pr.up[p1]
			B_value = pr.down[nX-p1]
		} else {
			p2a <- c(1:nX)[pr.down==tmax]
			if (length(p2a)>1) p2a <- p2a[1]

			A_value = pr.down[p2a]
			B_value = pr.up[nX-p2a]
		}
	} else if (tcmax==tmax) {
		ind<-row2
		col.part<-2
		if (ind==3) {
			A_value = pc.up[p3]
			B_value = pc.down[nY-p3]
		} else {
			p4a <- c(1:nY)[pc.down==tmax]
			if (length(p4a)>1) p4a <- p4a[1]
			
			A_value = pc.down[p4a]
			B_value = pc.up[nY-p4a]
		}
	} else {
		ind<-row3
		col.part<-3
		if (ind==5) {
			A_value = pt.up[p5]
			B_value = pt.down[nT-p5]
		} else {
			p6a <- c(1:nT)[pt.down==tmax]
			if (length(p6a)>1) p6a <- p6a[1]
			
			A_value = pt.down[p6a]
			B_value = pt.up[nT-p6a]
		}
	}
	
	# partition: [1] = initial index for row/col/tim partition
	#            [2] = last index for row/col/tim partition
	partition<-potential2[ind,]
	
	if (col.part==1) 
		gen.dim$nX<-partition
	else if (col.part==2) 
		gen.dim$nY<-partition 
	else
		gen.dim$nT<-partition 
	
	values <- list(A_value=A_value, B_value=B_value)
	
	# part: col.part = 1 if choosing row partition
	#                  2 if choosing col partition
	#                  3 if choosing tim partition
	#       partition = [1] = initial index for row/col/tim partition
	#                   [2] = last index for row/col/tim partition
	#       gen.dim: $nX & $nY, $nT updated with the chosen partition
	part <- list(col.part=col.part,partition=partition,gen.dim=gen.dim,values=values)
	
	part
}

####
# inputs, e.g. for rows:
#   cz = overall max z-score
#   pz = z-score for the last row index in upward direction
#   nsearch = total number of rows
#   cmu = expected sum for the rows with max z-score
#   pmu = overall expected sum for all rows and columns
#
# function: returning normalized z-score
#
partition.rule <- function(cz,pz,nsearch,cmu,pmu) {
  mu <- 0.0192357+pz*0.7842139+nsearch*0.0072527-nsearch^2*0.0001576-
  	cmu*0.0002191+pmu*0.0002004+
  	pz^2*0.4209625-pz*nsearch*0.0105873+pz*nsearch^2*0.0001593-
  	pz*cmu*0.0003409+pz*pmu*0.0003577-
  	cmu*pz^2*0.0003779+pmu*pz^2*0.0001439
  
  sigma <- exp(-3.3696115-pz*1.8448222+pz^2*0.7368984+nsearch*0.0333401-
  						 	cmu*0.0009042+pmu*0.0007915+
  						 	pz*nsearch*0.0574132-pz*cmu*0.0011366+pz*pmu*0.0012935-
  						 	pz^2*nsearch*0.0074938+pz^2*cmu*0.0003779+pz^2*pmu*0.0004892)
  
  Z <- (cz-mu)/sigma
  Z
}

####
# inputs:
#   col.part: = FALSE if choosing row partition
#               TRUE if choosing col partition
#   partition = [1] = initial index for row/col partition A
#               [2] = last index for row/col partition
#   scounts: the original matrix O to be analysed for row/col partitioning
#   smus: corresponding matrix with expected counts (same dimension as scounts)
#   gen.dim: $row: initial row and last row indices O2
#            $column: initial col and last col indices
#
# returning 2 offsprings A & B in ret:
#   hp.count = desired partition A, lp.count = other partition B
#   hp.mu, lp.mu: corresponding to hp.count, lp.count respectively
#   gen.dim1 = dimension A for hp.count, gen.dim2 = dimension B for lp.count
#
offspring<-function(col.part,partition,scount,smu,gen.dim,maxSum) {
  gen.dim1<-gen.dim
  gen.dim2<-gen.dim
  
  if(col.part==1) {
  	# hp.count = extracted desired row partition
  	# hp.mu = similar to hp.count of the same corresponding dimension
  	# gen.dim1 = extracted row, col & tim dimensions
  	hp.count <- scount[c(partition[1]:partition[2]),gen.dim$nY[1]:gen.dim$nY[2],
  										 gen.dim$nT[1]:gen.dim$nT[2]]
  	hp.mu <- smu[c(partition[1]:partition[2]),gen.dim$nY[1]:gen.dim$nY[2],
  							 gen.dim$nT[1]:gen.dim$nT[2]]
  	gen.dim1$nX <- partition
  	
  	# gen.dim2 = row & col dimensions of the other row partition
  	if (partition[1]==gen.dim$nX[1]) {
  		if (partition[2]!=gen.dim$nX[2])
  			gen.dim2$nX<-c(partition[2]+1,gen.dim$nX[2])
  		else
  			# should not arrive here anyway
  			gen.dim2$nX<-c(gen.dim$nX[2],gen.dim$nX[2])
  	} else
  		gen.dim2$nX<-c(gen.dim$nX[1],partition[1]-1)
  } else if (col.part==2) {
    # hp.count = extracted desired col partition
    # hp.mu = similar to hp.count of the same corresponding dimension
    # gen.dim1 = extracted row, col & tim dimensions
    hp.count <- scount[gen.dim$nX[1]:gen.dim$nX[2],partition[1]:partition[2],
    									 gen.dim$nT[1]:gen.dim$nT[2]]
    hp.mu <- smu[gen.dim$nX[1]:gen.dim$nX[2],partition[1]:partition[2],
    						 gen.dim$nT[1]:gen.dim$nT[2]]
    gen.dim1$nY<-partition
    
    # gen.dim2 = row & col dimensions of the other col partition
    if (partition[1]==gen.dim$nY[1]) {
      if (partition[2]!=gen.dim$nY[2])
      	gen.dim2$nY<-c(partition[2]+1,gen.dim$nY[2])
      else
      	# should not arrive here anyway
      	gen.dim2$nY<-c(gen.dim$nY[2],gen.dim$nY[2])
    } else
    	gen.dim2$nY<-c(gen.dim$nY[1],partition[1]-1)
  } else {
  	# hp.count = extracted desired tim partition
  	# hp.mu = similar to hp.count of the same corresponding dimension
  	# gen.dim1 = extracted row, col & tim dimensions
    hp.count <- scount[gen.dim$nX[1]:gen.dim$nX[2],gen.dim$nY[1]:gen.dim$nY[2],
    									 c(partition[1]:partition[2])];
    hp.mu <- smu[gen.dim$nX[1]:gen.dim$nX[2],gen.dim$nY[1]:gen.dim$nY[2],
    						 c(partition[1]:partition[2])]
    gen.dim1$nT <- partition
    
  	# gen.dim2 = row & col dimensions of the other row partition
  	if (partition[1]==gen.dim$nT[1]) {
      if (partition[2]!=gen.dim$nT[2])
      	gen.dim2$nT<-c(partition[2]+1,gen.dim$nT[2])
      else
      	# should not arrive here anyway
      	gen.dim2$nT<-c(gen.dim$nT[2],gen.dim$nT[2])
    } else
    	gen.dim2$nT<-c(gen.dim$nT[1],partition[1]-1)
  }
  # lp.count = the other partition
  # lp.mu = similar to lp.count of the same corresponding dimension
  lp.count<-scount[gen.dim2$nX[1]:gen.dim2$nX[2],gen.dim2$nY[1]:gen.dim2$nY[2],
  								 gen.dim2$nT[1]:gen.dim2$nT[2]]
  lp.mu<-smu[gen.dim2$nX[1]:gen.dim2$nX[2],gen.dim2$nY[1]:gen.dim2$nY[2],
  					 gen.dim2$nT[1]:gen.dim2$nT[2]]
  
  if (sum(hp.count) <= maxSum) StopA = TRUE
  else StopA = FALSE
  
  if (sum(lp.count) <= maxSum) StopB = TRUE
  else StopB = FALSE
  
  # hp.count = desired partition, lp.count = other partition
  # hp.mu, lp.mu: corresponding to hp.count, lp.count respectively
  # gen.dim1 = dimension for hp.count, gen.dim2 = dimension for lp.count
  # STOP = TRUE if desired partition area same as original area
  ret <- list(hp.count=hp.count,lp.count=lp.count,
  						hp.mu=hp.mu,lp.mu=lp.mu,
  						gen.dim1=gen.dim1,gen.dim2=gen.dim2,StopA=StopA,StopB=StopB)
  
  ret
}

copyParentA <- function(gen.dim) {
	gen.dim1 <- gen.dim
	gen.dim2 <- list(nX=c(0,0),nY=c(0,0),nT=c(0,0))
	
	hp.count <- NULL; hp.mu <- NULL
	lp.count <- NULL; lp.mu <- NULL

	StopA = TRUE; StopB = TRUE
	
	ret <- list(hp.count=hp.count,lp.count=lp.count,
							hp.mu=hp.mu,lp.mu=lp.mu,
							gen.dim1=gen.dim1,gen.dim2=gen.dim2,StopA=StopA,StopB=StopB)
	ret
}

copyParentB <- function(gen.dim) {
	gen.dim1 <- list(nX=c(0,0),nY=c(0,0),nT=c(0,0))
	gen.dim2 <- gen.dim
	
	hp.count <- NULL; hp.mu <- NULL
	lp.count <- NULL; lp.mu <- NULL
	
	StopA = TRUE; StopB = TRUE
	
	ret <- list(hp.count=hp.count,lp.count=lp.count,
							hp.mu=hp.mu,lp.mu=lp.mu,
							gen.dim1=gen.dim1,gen.dim2=gen.dim2,StopA=StopA,StopB=StopB)
	ret
}

####

# inputs:
#       scounts: the original matrix O to be analysed for row/col partitioning
#       smus: corresponding matrix with expected counts (same dimension as scounts)
#       gen1.dim1: dimension for desired partition A
#       gen1.dim2: dimension for other partition B
#       num: number of iterations (for reference, unused yet)
#
# returning 4 offsprings A-A, A-B, B-A, B-B in tem:
#   hp.g21.count = partition A-A, lp.g21.count = partition A-B
#   hp.g21.mu lp.g21.mu: expected counts for partitions A-A, A-B respectively
#   hp.g22.count = partition B-A, lp.g22.count = partition B-B
#   hp.g22.mu lp.g22.mu: expected counts for partitions B-A, B-B respectively
#   gen21.dim1, gen21.dim2, gen22.dim1, gen22.dim2 = dimensions for A-A, A-B, B-A, B-B
#   STOP1 = TRUE if A-A same as A, STOP2 = TRUE if B-A same as B
next.generation <- function(scount,smu,gen1.dim1,gen1.dim2,num) {
	# choosing the best row/col partition A-A for 1 iteration:
  part <- recursive.partition(
  	scount[gen1.dim1$row[1]:gen1.dim1$row[2],gen1.dim1$column[1]:gen1.dim1$column[2]],
  	   smu[gen1.dim1$row[1]:gen1.dim1$row[2],gen1.dim1$column[1]:gen1.dim1$column[2]],
  	   gen1.dim1, num)
  
  #scount<-scount
  #smu<-smu
  
  # gen1.dim.ng = dimension for partition A
  gen1.dim.ng<-gen1.dim1
  
  if (part$partition[1]>part$partition[2]) STOP7
  if (part$col.part==TRUE && part$partition[2]>gen1.dim1$column[2]) STOP8
  
  # numm: for reference, unused yet
  numm <- num+1
  
  # children.gen21 = get 2 offsprings A-A and A-B
  children.gen21<-offspring(part$col.part,part$partition,scount,smu,gen1.dim1,num)
  
  # numm: for reference, unused yet
  numm2<-numm+1
  
  if(children.gen21$gen.dim1$column[1]>children.gen21$gen.dim1$column[2])STOP4
  if(children.gen21$gen.dim2$column[1]>children.gen21$gen.dim2$column[2])STOP6
  
  # choosing the best row/col partition B-A for 1 iteration:
  part <- recursive.partition(
  	scount[gen1.dim2$row[1]:gen1.dim2$row[2],gen1.dim2$column[1]:gen1.dim2$column[2]],
       smu[gen1.dim2$row[1]:gen1.dim2$row[2],gen1.dim2$column[1]:gen1.dim2$column[2]],
  	   gen1.dim2,num)
  
  if(part$partition[1]>part$partition[2])STOP9
  if(part$col.part==TRUE && part$partition[2]>gen1.dim2$column[2])STOP10
  
  # children.gen22 = get 2 offsprings B-A and B-B
  children.gen22<-offspring(part$col.part,part$partition,scount,smu,gen1.dim2,num)
  
  # numm3: for reference, unused yet
  numm3<-numm2+1
  
  if(children.gen22$gen.dim1$column[1]>children.gen22$gen.dim1$column[2])STOP5
  if(children.gen22$gen.dim2$column[1]>children.gen22$gen.dim2$column[2])STOP7
  
  # hp.g21.count = partition A-A, lp.g21.count = partition A-B
  # hp.g21.mu lp.g21.mu: expected counts for partitions A-A, A-B respectively
  # hp.g22.count = partition B-A, lp.g22.count = partition B-B
  # hp.g22.mu lp.g22.mu: expected counts for partitions B-A, B-B respectively
  # gen21.dim1, gen21.dim2, gen22.dim1, gen22.dim2 = dimensions for A-A, A-B, B-A, B-B
  # STOP1 = TRUE if A-A same as A, STOP2 = TRUE if B-A same as B
  tem<-list(hp.g21.count=children.gen21$hp.count, lp.g21.count=children.gen21$lp.count,
  					hp.g21.mu=children.gen21$hp.mu, lp.g21.mu=children.gen21$lp.mu,
            hp.g22.count=children.gen22$hp.count, lp.g22.count=children.gen22$lp.count,
  					hp.g22.mu=children.gen22$hp.mu,lp.g22.mu=children.gen22$lp.mu,
            gen21.dim1=children.gen21$gen.dim1, gen21.dim2=children.gen21$gen.dim2,
  					gen22.dim1=children.gen22$gen.dim1, gen22.dim2=children.gen22$gen.dim2,
  					STOP1=children.gen21$STOP,STOP2=children.gen22$STOP)
  
  tem<-tem
  tem
}

#### 
# inputs:
#   counts: 2D observed counts
#   mu: 2D expected means
#
# returning zmax = the overall max z score after all row & col partitioning
#
group.rp<-function(counts,mu) {
  tempmu<-mu
  temp<-counts

  N = length(tempmu[1,])
  
  zmax<-2*(sqrt(sum(temp))-sqrt(sum(tempmu)))
  
  gen.dim<-list(row=c(1,N),column=c(1,N))
  
  gen1.dim1<-gen.dim
  gen1.dim2<-gen.dim
  
  part<-recursive.partition(temp,tempmu,gen.dim,1)
  
  if ((part$partition[1]-1+N-part$partition[2])>0) {
    children.gen1 <- offspring(part$col.part,part$partition,temp,tempmu,gen.dim,1)
    
    if (children.gen1$STOP==FALSE) {
      children.gen2 <- next.generation(temp,tempmu,
      									 children.gen1$gen.dim1,children.gen1$gen.dim2,1)
      zmax <- max(zmax,2*(sqrt(sum(children.gen1$hp.count))-sqrt(sum(children.gen1$hp.mu))))
      
      if(children.gen2$STOP1==FALSE) {
        children.gen31 <- next.generation(temp,tempmu,
        										children.gen2$gen21.dim1,children.gen2$gen21.dim2,2)
        zmax <- max(zmax,c(
        	2*(sqrt(sum(children.gen31$hp.g21.count))-sqrt(sum(children.gen31$hp.g21.mu))),
        	2*(sqrt(sum(children.gen31$hp.g22.count))-sqrt(sum(children.gen31$hp.g22.mu)))))
        
        if(children.gen31$STOP1==FALSE) {
          children.gen411 <- next.generation(temp,tempmu,
          										 children.gen31$gen21.dim1,children.gen31$gen21.dim2,3)
          zmax<-max(zmax,c(
          	2*(sqrt(sum(children.gen411$hp.g21.count))-sqrt(sum(children.gen411$hp.g21.mu))),
          	2*(sqrt(sum(children.gen411$hp.g22.count))-sqrt(sum(children.gen411$hp.g22.mu)))))
          
          if(children.gen411$STOP1==FALSE) {
            children.gen5111 <- next.generation(temp,tempmu,
            			                children.gen411$gen21.dim1,children.gen411$gen21.dim2,4)
            zmax <- max(zmax,c(
            	2*(sqrt(sum(children.gen5111$hp.g21.count))-sqrt(sum(children.gen5111$hp.g21.mu))),
            	2*(sqrt(sum(children.gen5111$hp.g22.count))-sqrt(sum(children.gen5111$hp.g22.mu)))))
            
            if(children.gen5111$STOP1==FALSE) {
              children.gen61111 <- next.generation(temp,tempmu,
              											 children.gen5111$gen21.dim1,children.gen5111$gen21.dim2,5)
              zmax<-max(zmax,c(
              	2*(sqrt(sum(children.gen61111$hp.g21.count))-
              		 	sqrt(sum(children.gen61111$hp.g21.mu))),
              	2*(sqrt(sum(children.gen61111$hp.g22.count))-
              		 	sqrt(sum(children.gen61111$hp.g22.mu)))))
            }
            
            if(children.gen5111$STOP2==FALSE) {
              children.gen61112 <- next.generation(temp,tempmu,
              											 children.gen5111$gen22.dim1,children.gen5111$gen22.dim2,6)
              zmax <- max(zmax,
              	2*(sqrt(sum(children.gen61112$hp.g21.count))-
              		 	sqrt(sum(children.gen61112$hp.g21.mu))),
              	2*(sqrt(sum(children.gen61112$hp.g22.count))-
              		 	sqrt(sum(children.gen61112$hp.g22.mu))))
            }
          }
          
          if(children.gen411$STOP2==FALSE){
            children.gen5112 <- next.generation(temp,tempmu,
            											children.gen411$gen22.dim1,children.gen411$gen22.dim2,7)
            zmax<-max(zmax,
            	2*(sqrt(sum(children.gen5112$hp.g21.count))-sqrt(sum(children.gen5112$hp.g21.mu))),
            	2*(sqrt(sum(children.gen5112$hp.g22.count))-sqrt(sum(children.gen5112$hp.g22.mu))))
            
            if (children.gen5112$STOP1==FALSE){
              children.gen61121 <- next.generation(temp,tempmu,
              											 children.gen5112$gen21.dim1,children.gen5112$gen21.dim2,8)
              zmax <- max(zmax,
              	2*(sqrt(sum(children.gen61121$hp.g21.count))-
              		 	sqrt(sum(children.gen61121$hp.g21.mu))),
                2*(sqrt(sum(children.gen61121$hp.g22.count))-
                	 	sqrt(sum(children.gen61121$hp.g22.mu))))
            }
            
            if (children.gen5112$STOP2==FALSE){
              children.gen61122 <- next.generation(temp,tempmu,
              											 children.gen5112$gen22.dim1,children.gen5112$gen22.dim2,9)
              zmax <- max(zmax,
              	2*(sqrt(sum(children.gen61122$hp.g21.count))-
              	    sqrt(sum(children.gen61122$hp.g21.mu))),
              	2*(sqrt(sum(children.gen61122$hp.g22.count))-
              		 	sqrt(sum(children.gen61122$hp.g22.mu))))}
          }
        }
        
        if (children.gen31$STOP2==FALSE) {
          children.gen412 <- next.generation(temp,tempmu,
          										 children.gen31$gen22.dim1,children.gen31$gen22.dim2,10)
          zmax<-max(zmax,
          	2*(sqrt(sum(children.gen412$hp.g21.count))-
          		  sqrt(sum(children.gen412$hp.g21.mu))),
          	2*(sqrt(sum(children.gen412$hp.g22.count))-
          		 	sqrt(sum(children.gen412$hp.g22.mu))))
          
          if (children.gen412$STOP1==FALSE) {
            children.gen5121<-next.generation(temp,tempmu,
            										children.gen412$gen21.dim1,children.gen412$gen21.dim2,11)
            zmax <- max(zmax,
            	2*(sqrt(sum(children.gen5121$hp.g21.count))-
            	    sqrt(sum(children.gen5121$hp.g21.mu))),
            	2*(sqrt(sum(children.gen5121$hp.g22.count))-
            		  sqrt(sum(children.gen5121$hp.g22.mu))))
            
            if(children.gen5121$STOP1==FALSE){
              children.gen61211 <- next.generation(temp,tempmu,
              										   children.gen5121$gen21.dim1,children.gen5121$gen21.dim2,12)
              zmax<-max(zmax,
              	2*(sqrt(sum(children.gen61211$hp.g21.count))-
              		  sqrt(sum(children.gen61211$hp.g21.mu))),
              	2*(sqrt(sum(children.gen61211$hp.g22.count))-
              		  sqrt(sum(children.gen61211$hp.g22.mu))))
            }
            
            if(children.gen5121$STOP2==FALSE){
              children.gen61212 <- next.generation(temp,tempmu,
              											 children.gen5121$gen22.dim1,children.gen5121$gen22.dim2,13)
              zmax<-max(zmax,
              	2*(sqrt(sum(children.gen61212$hp.g21.count))-
              		 	sqrt(sum(children.gen61212$hp.g21.mu))),
              	2*(sqrt(sum(children.gen61212$hp.g22.count))-
              		 	sqrt(sum(children.gen61212$hp.g22.mu))))
            }  
          }
          
          if (children.gen412$STOP2==FALSE) {
            children.gen5122<-next.generation(temp,tempmu,
            										children.gen412$gen22.dim1,children.gen412$gen22.dim2,14) 
            zmax<-max(zmax,
            	2*(sqrt(sum(children.gen5122$hp.g21.count))-
            		 	sqrt(sum(children.gen5122$hp.g21.mu))),
            	2*(sqrt(sum(children.gen5122$hp.g22.count))-
            		 	sqrt(sum(children.gen5122$hp.g22.mu))))
            
            if (children.gen5122$STOP1==FALSE) {
              children.gen61221<-next.generation(temp,tempmu,
              										 children.gen5122$gen21.dim1,children.gen5122$gen21.dim2,15)
              zmax<-max(zmax,
              	2*(sqrt(sum(children.gen61221$hp.g21.count))-
              		 	sqrt(sum(children.gen61221$hp.g21.mu))),
              	2*(sqrt(sum(children.gen61221$hp.g22.count))-
              		 	sqrt(sum(children.gen61221$hp.g22.mu))))
            }
            
            if (children.gen5122$STOP2==FALSE) {
              children.gen61222<-next.generation(temp,tempmu,
              										 children.gen5122$gen22.dim1,children.gen5122$gen22.dim2,16)
              zmax<-max(zmax,
              	2*(sqrt(sum(children.gen61222$hp.g21.count))-
              		 	sqrt(sum(children.gen61222$hp.g21.mu))),
              	2*(sqrt(sum(children.gen61222$hp.g22.count))-
              		 	sqrt(sum(children.gen61222$hp.g22.mu))))
            }
          }
        }
      }
      
      if(children.gen2$STOP2==FALSE) {
        children.gen32 <- next.generation(temp,tempmu,
        										children.gen2$gen22.dim1,children.gen2$gen22.dim2,17)
        zmax<-max(zmax,
        	2*(sqrt(sum(children.gen32$hp.g21.count))-
        		 	sqrt(sum(children.gen32$hp.g21.mu))),
        	2*(sqrt(sum(children.gen32$hp.g22.count))-
        		 	sqrt(sum(children.gen32$hp.g22.mu))))
        
        if(children.gen32$STOP1==FALSE) {
          children.gen421<-next.generation(temp,tempmu,
          									 children.gen32$gen21.dim1,children.gen32$gen21.dim2,18)
          zmax<-max(zmax,
          	2*(sqrt(sum(children.gen421$hp.g21.count))-
          		 	sqrt(sum(children.gen421$hp.g21.mu))),
          	2*(sqrt(sum(children.gen421$hp.g22.count))-
          		 	sqrt(sum(children.gen421$hp.g22.mu))))
          
          if(children.gen421$STOP1==FALSE) {
            children.gen5211<-next.generation(temp,tempmu,
            										children.gen421$gen21.dim1,children.gen421$gen21.dim2,19)
            zmax<-max(zmax,
            	2*(sqrt(sum(children.gen5211$hp.g21.count))-
            		 	sqrt(sum(children.gen5211$hp.g21.mu))),
            	2*(sqrt(sum(children.gen5211$hp.g22.count))-
            		 	sqrt(sum(children.gen5211$hp.g22.mu))))
            
            if(!is.null(children.gen5211$STOP1)) {
            	if(children.gen5211$STOP1==FALSE) {
                children.gen62111<-next.generation(temp,tempmu,
                										 children.gen5211$gen21.dim1,children.gen5211$gen21.dim2,20)
	              zmax<-max(zmax,
	              	2*(sqrt(sum(children.gen62111$hp.g21.count))-
	              		 	sqrt(sum(children.gen62111$hp.g21.mu))),
	              	2*(sqrt(sum(children.gen62111$hp.g22.count))-
	              		 	sqrt(sum(children.gen62111$hp.g22.mu))))}
            }
            
            if (!is.null(children.gen5211$STOP2)) { 
            	if (!is.null(children.gen5211$STOP2)) { 
            		if (children.gen5211$STOP2==FALSE) {
              		children.gen62112<-next.generation(temp,tempmu,
              											 children.gen5211$gen22.dim1,children.gen5211$gen22.dim2,21)
		              zmax<-max(zmax,
		              	2*(sqrt(sum(children.gen62112$hp.g21.count))-
		              		 	sqrt(sum(children.gen62112$hp.g21.mu))),
		              	2*(sqrt(sum(children.gen62112$hp.g22.count))-
		              		 	sqrt(sum(children.gen62112$hp.g22.mu))))
            		}
            	}
            }
          }
          
          if(children.gen421$STOP2==FALSE) {
            children.gen5212<-next.generation(temp,tempmu,
            										children.gen421$gen22.dim1,children.gen421$gen22.dim2,22)
            zmax<-max(zmax,
            	2*(sqrt(sum(children.gen5212$hp.g21.count))-
            		 	sqrt(sum(children.gen5212$hp.g21.mu))),
            	2*(sqrt(sum(children.gen5212$hp.g22.count))-
            		 	sqrt(sum(children.gen5212$hp.g22.mu))))
            
            if(!is.null(children.gen5212$STOP1)) { 
            	if(children.gen5212$STOP1==FALSE) {
                children.gen62121<-next.generation(temp,tempmu,
                									 children.gen5212$gen21.dim1,children.gen5212$gen21.dim2,23)
	              zmax<-max(zmax,
	              	2*(sqrt(sum(children.gen62121$hp.g21.count))-
	              		 	sqrt(sum(children.gen62121$hp.g21.mu))),
	              	2*(sqrt(sum(children.gen62121$hp.g22.count))-
	              		 	sqrt(sum(children.gen62121$hp.g22.mu))))
            	}
            }
            
            if(!is.null(children.gen5212$STOP2)) { 
            	if(children.gen5212$STOP2==FALSE) {
                children.gen62122<-next.generation(temp,tempmu,
                									 children.gen5212$gen22.dim1,children.gen5212$gen22.dim2,24)
              zmax<-max(zmax,
              	2*(sqrt(sum(children.gen62122$hp.g21.count))-
              		 	sqrt(sum(children.gen62122$hp.g21.mu))),
              	2*(sqrt(sum(children.gen62122$hp.g22.count))-
              		 	sqrt(sum(children.gen62122$hp.g22.mu))))
            	}
            }
          } 
        }
        
        if (!is.null(children.gen32$STOP1)) { 
        	if (children.gen32$STOP2==FALSE) {
            children.gen422<-next.generation(temp,tempmu,
            								 children.gen32$gen22.dim1,children.gen32$gen22.dim2,25)
	          zmax<-max(zmax,
	          	2*(sqrt(sum(children.gen422$hp.g21.count))-
	          		 	sqrt(sum(children.gen422$hp.g21.mu))),
	          	2*(sqrt(sum(children.gen422$hp.g22.count))-
	          		 	sqrt(sum(children.gen422$hp.g22.mu))))
	          
	          if(children.gen422$STOP1==FALSE) {
	            children.gen5221<-next.generation(temp,tempmu,
	            							    children.gen422$gen21.dim1,children.gen422$gen21.dim2,26)
	            zmax<-max(zmax,
	            	2*(sqrt(sum(children.gen5221$hp.g21.count))-
	            		 	sqrt(sum(children.gen5221$hp.g21.mu))),
	            	2*(sqrt(sum(children.gen5221$hp.g22.count))-
	            		 	sqrt(sum(children.gen5221$hp.g22.mu))))
	            
	            if(children.gen5221$STOP1==FALSE) {
	              children.gen62211<-next.generation(temp,tempmu,
	              									 children.gen5221$gen21.dim1,children.gen5221$gen21.dim2,27)
	              zmax<-max(zmax,
	              	2*(sqrt(sum(children.gen62211$hp.g21.count))-
	              		 	sqrt(sum(children.gen62211$hp.g21.mu))),
	              	2*(sqrt(sum(children.gen62211$hp.g22.count))-
	              		 	sqrt(sum(children.gen62211$hp.g22.mu))))
	            }
	            
	            if(children.gen5221$STOP2==FALSE){
	              children.gen62212<-next.generation(temp,tempmu,
	              									 children.gen5221$gen22.dim1,children.gen5221$gen22.dim2,28)
	              zmax<-max(zmax,
	              	2*(sqrt(sum(children.gen62212$hp.g21.count))-
	              		 	sqrt(sum(children.gen62212$hp.g21.mu))),
	              	2*(sqrt(sum(children.gen62212$hp.g22.count))-
	              		 	sqrt(sum(children.gen62212$hp.g22.mu))))
	            }
	          }
	          
	          if(children.gen422$STOP2==FALSE){
	            children.gen5222<-next.generation(temp,tempmu,
	            									children.gen422$gen22.dim1,children.gen422$gen22.dim2,29)
	            zmax<-max(zmax,
	            	2*(sqrt(sum(children.gen5222$hp.g21.count))-
	            		 	sqrt(sum(children.gen5222$hp.g21.mu))),
	            	2*(sqrt(sum(children.gen5222$hp.g22.count))-
	            		 	sqrt(sum(children.gen5222$hp.g22.mu))))
	            
	            if(children.gen5222$STOP1==FALSE) {
	              children.gen62221<-next.generation(temp,tempmu,
	              									 children.gen5222$gen21.dim1,children.gen5222$gen21.dim2,30)
	              zmax<-max(zmax,
	              	2*(sqrt(sum(children.gen62221$hp.g21.count))-
	              		 	sqrt(sum(children.gen62221$hp.g21.mu))),
	              	2*(sqrt(sum(children.gen62221$hp.g22.count))-
	              		 	sqrt(sum(children.gen62221$hp.g22.mu))))
	            }
	            
	            if(children.gen5222$STOP1==FALSE) {
	              children.gen62222<-next.generation(temp,tempmu,
	              									 children.gen5222$gen22.dim1,children.gen5222$gen22.dim2,31)
	              zmax<-max(zmax,
	              	2*(sqrt(sum(children.gen62222$hp.g22.count))-
	              		 	sqrt(sum(children.gen62222$hp.g22.mu))),
	              	2*(sqrt(sum(children.gen62222$hp.g21.count))-
	              		 	sqrt(sum(children.gen62222$hp.g21.mu))))
	            }
	          }
	        }
	      }
  	  }
  	}
  }
  
  zmax
}

####

appendList <- function (x, val) 
{
  stopifnot(is.list(x), is.list(val))
  xnames <- names(x)
  
  for (v in names(val)) {
    x[[v]] <- 
    	if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]])) 
        appendList(x[[v]], val[[v]])
      else 
      	c(x[[v]], val[[v]])
  }
  
  x
}

####

recordData <- function(n, c, f, children, p)
{
	dim1 = children$gen.dim1
	dim2 = children$gen.dim2
	
	myData <- data.frame(iteration=n, child=c, col.part=f, 
											 xA1=dim1$nX[1], xA2=dim1$nX[2], 
											 yA1=dim1$nY[1], yA2=dim1$nY[2], 
											 tA1=dim1$nT[1], tA2=dim1$nT[2], 
											 StopA=children$StopA, A_value=p$A_value,
											 xB1=dim2$nX[1], xB2=dim2$nX[2],
											 yB1=dim2$nY[1], yB2=dim2$nY[2], 
											 tB1=dim2$nT[1], tB2=dim2$nT[2], 
											 StopB=children$StopB, B_value=p$B_value)
	return (myData)
}

####

getData <- function(counts, mu, curData1, index)
{
	A.dim <- NULL
	A.counts <- NULL
	A.mu <- NULL
	A.Stop <- FALSE

	B.dim <- NULL
	B.counts <- NULL
	B.mu <- NULL
	B.Stop <- FALSE
	
	A_value <- NULL
	B_value <- NULL
	
	maxRow = dim(curData1)[1]
	if (maxRow >= index)
	{
		row1 = curData1$xA1[index]
		row2 = curData1$xA2[index]
		col1 = curData1$yA1[index]
		col2 = curData1$yA2[index]
		tim1 = curData1$tA1[index]
		tim2 = curData1$tA2[index]
		
		A.dim <- list(nX=c(row1,row2),nY=c(col1,col2),nT=c(tim1,tim2))

		dim = c(row2-row1+1,col2-col1+1,tim2-tim1+1)
		A.counts <- array(counts[row1:row2, col1:col2, tim1:tim2], dim)
		A.mu <- array(mu[row1:row2, col1:col2, tim1:tim2], dim)
		
		A.Stop <- curData1$StopA[index]
		if (prod(dim) == 1) A.Stop <- TRUE
		
		A_value <- curData1$A_value[index]

		row1 = curData1$xB1[index]
		row2 = curData1$xB2[index]
		col1 = curData1$yB1[index]
		col2 = curData1$yB2[index]
		tim1 = curData1$tB1[index]
		tim2 = curData1$tB2[index]
		
		B.dim <- list(nX=c(row1,row2),nY=c(col1,col2),nT=c(tim1,tim2))
		
		dim = c(row2-row1+1,col2-col1+1,tim2-tim1+1)
		B.counts <- array(counts[row1:row2, col1:col2, tim1:tim2], dim)
		B.mu <- array(mu[row1:row2, col1:col2, tim1:tim2], dim)
		
		B.Stop <- curData1$StopB[index]
		if (prod(dim) == 1) B.Stop <- TRUE
		
		B_value <- curData1$B_value[index]
	}
	
	retData <- list(A.counts = A.counts, A.mu = A.mu, A.dim=A.dim, 
									A.Stop = A.Stop, A_value=A_value,
									B.counts = B.counts, B.mu = B.mu, B.dim=B.dim, 
									B.Stop = B.Stop, B_value=B_value)
	
	return (retData)
}

formatNum <- function(n)
{
	if (n == "")
		num1 <- ""
	else
	{
		num1 <- as.numeric(n)
		num1 <- format(num1, digits=3)
	}
	return (num1)
}

getValues <- function(m, ite)
{
	# nrow = dim(myData)[1]

	if (ite==1)
	{
		curData2 <- data.frame(xA1=m$xA1[1], xA2=m$xA2[1], 
													 yA1=m$yA1[1], yA2=m$yA2[1],
													 tA1=m$tA1[1], tA2=m$tA2[1], 
													 A_value=formatNum(m$A_value[1]), StopA=m$StopA[1],
													 xB1=m$xB1[1], xB2=m$xB2[1], 
													 yB1=m$yB1[1], yB2=m$yB2[1],
													 tB1=m$tB1[1], tB2=m$tB2[1], 
													 B_value=formatNum(m$B_value[1]), StopB=m$StopB[1])
	} else {
		index <- 0
		for (i in 1:(ite-1))
		{
			index <- index + 2^(i-1)
		}

		curData2 <- data.frame(xA1=m$xA1[index+1], xA2=m$xA2[index+1], 
													 yA1=m$yA1[index+1], yA2=m$yA2[index+1],
													 tA1=m$tA1[index+1], tA2=m$tA2[index+1], 
													 A_value=formatNum(m$A_value[index+1]), StopA=m$StopA[index+1],
													 xB1=m$xB1[index+1], xB2=m$xB2[index+1], 
													 yB1=m$yB1[index+1], yB2=m$yB2[index+1],
													 tB1=m$tB1[index+1], tB2=m$tB2[index+1], 
													 B_value=formatNum(m$B_value[index+1]), StopB=m$StopB[index+1])
		for (child in 2:2^(ite-1))
		{
			i <- index + child
			tempData <- data.frame(xA1=m$xA1[i], xA2=m$xA2[i], 
														 yA1=m$yA1[i], yA2=m$yA2[i],
														 tA1=m$tA1[i], tA2=m$tA2[i], 
														 A_value=formatNum(m$A_value[i]), StopA=m$StopA[i],
														 xB1=m$xB1[i], xB2=m$xB2[i], 
														 yB1=m$yB1[i], yB2=m$yB2[i],
														 tB1=m$tB1[i], tB2=m$tB2[i], 
														 B_value=formatNum(m$B_value[i]), StopB=m$StopB[i])
			curData2 <- rbind(curData2, tempData)
		}
	}

	return (curData2)
}

plotA <- function(x, y, p, index, colA, colP, displayT)
{
	if (displayT>0)
	{
		if ((displayT < p$tA1[index]) || (p$tA2[index] < displayT)) return ()
	}
	if (p$xA1[index] != 0) {
		colorA <- ifelse(p$StopA[index],colP,colA)
		text(x,y,p$A_value[index],cex=0.8,col=colorA)
	}
}

plotA2 <- function(x, y, p, index, colA, colP, displayT)
{
	if (displayT>0)
	{
		if ((displayT < p$tA1[index]) || (p$tA2[index] < displayT)) return ()
	}
	if (p$xA1[index] != 0) {
		colorA <- ifelse(p$StopA[index],colP,colA)
		text(x,y,p$A_value[index],cex=0.8,col=colorA,pos=4)
	}
}

plotAWithValues <- function(x, y, p, index, colA, colP, displayT)
{
	if (y<0) return (FALSE)
	if (displayT>0)
	{
		if ((displayT < p$tA1[index]) || (p$tA2[index] < displayT)) return (FALSE)
	}
	if (p$StopA[index]) return (FALSE)
	if (p$xA1[index] != 0) {
		colorA <- ifelse(p$StopA[index],colP,colA)
		str <- paste(p$A_value[index],"  (", 
								 as.character(p$xA1[index]), ":", as.character(p$xA2[index]), ",",
								 as.character(p$yA1[index]), ":", as.character(p$yA2[index]), ")  [t=",
								 as.character(p$tA1[index]), ":", as.character(p$tA2[index]), "]", sep=" ")
		text(x,y,str,cex=0.8,col=colorA,pos=4)
		
		return (TRUE)
	} else
		return (FALSE)
}

plotB <- function(x, y, p, index, colB, colP, displayT)
{
	if (displayT>0)
	{
		if ((displayT < p$tB1[index]) || (p$tB2[index] < displayT)) return ()
	}
	if (p$xB1[index] != 0) {
		colorB <- ifelse(p$StopB[index],colP,colB)
		text(x,y,p$B_value[index],cex=0.8,col=colorB)
	}
}

plotB2 <- function(x, y, p, index, colB, colP, displayT)
{
	if (displayT>0)
	{
		if ((displayT < p$tB1[index]) || (p$tB2[index] < displayT)) return ()
	}
	if (p$xB1[index] != 0) {
		colorB <- ifelse(p$StopB[index],colP,colB)
		text(x,y,p$B_value[index],cex=0.8,col=colorB,pos=4)
	}
}

plotBWithValues <- function(x, y, p, index, colB, colP, displayT)
{
	if (y<0) return (FALSE)
	if (displayT>0)
	{
		if ((displayT < p$tB1[index]) || (p$tB2[index] < displayT)) return (FALSE)
	}
	if (p$StopB[index]) return (FALSE)
	if (p$xB1[index] != 0) {
		colorB <- ifelse(p$StopB[index],colP,colB)
		str <- paste(p$B_value[index],"  (", 
								 as.character(p$xB1[index]), ":", as.character(p$xB2[index]), ",",
								 as.character(p$yB1[index]), ":", as.character(p$yB2[index]), ")  [t=",
								 as.character(p$tB1[index]), ":", as.character(p$tB2[index]), "]", sep=" ")
		text(x,y,str,cex=0.8,col=colorB,pos=4)
		
		return (TRUE)
	} else
		return (FALSE)
}

plot1Row <- function(y, p, index, colA, colB, colP, displayT)
{
	plotA(1,y,p,index,colA,colP,displayT)
	plotB(4,y,p,index,colB,colP,displayT)
	plotA(7,y,p,index+1,colA,colP,displayT)
	plotB(10,y,p,index+1,colB,colP,displayT)
}

display_values <- function(p, ite, colA, colB, colP, nrow, displayT)
{
	h <- nrow+1   # height
	
	if (ite==1) {
		y <- h*0.5
		plotA(1+1.5,y,p,1,colA,colP,displayT)
		plotB(7+1.5,y,p,1,colB,colP,displayT)
	}
	else if (ite==2) {
		y <- h*0.5; plot1Row(y,p,1,colA,colB,colP,displayT)
	}
	else if (ite==3) {
		y <- h*0.75; plot1Row(y,p,1,colA,colB,colP,displayT)
		y <- h*0.25; plot1Row(y,p,3,colA,colB,colP,displayT)
	}
	else if (ite==4) {
		y <- h*0.8; plot1Row(y,p,1,colA,colB,colP,displayT)
		y <- h*0.6; plot1Row(y,p,3,colA,colB,colP,displayT)
		y <- h*0.4; plot1Row(y,p,5,colA,colB,colP,displayT)
		y <- h*0.2; plot1Row(y,p,7,colA,colB,colP,displayT)
	}
	else if (ite==5) {
		y <- h*0.8; plot1Row(y,p,1,colA,colB,colP,displayT)
		y <- h*0.7; plot1Row(y,p,3,colA,colB,colP,displayT)
		y <- h*0.6; plot1Row(y,p,5,colA,colB,colP,displayT)
		y <- h*0.5; plot1Row(y,p,7,colA,colB,colP,displayT)

		y <- h*0.4; plot1Row(y,p,9,colA,colB,colP,displayT)
		y <- h*0.3; plot1Row(y,p,11,colA,colB,colP,displayT)
		y <- h*0.2; plot1Row(y,p,13,colA,colB,colP,displayT)
		y <- h*0.1; plot1Row(y,p,15,colA,colB,colP,displayT)
	}
	else if (ite==6) {
		y <- h*15/16; plot1Row(y,p,1,colA,colB,colP,displayT)
		y <- h*14/16; plot1Row(y,p,3,colA,colB,colP,displayT)
		y <- h*13/16; plot1Row(y,p,5,colA,colB,colP,displayT)
		y <- h*12/16; plot1Row(y,p,7,colA,colB,colP,displayT)
		
		y <- h*11/16; plot1Row(y,p,9,colA,colB,colP,displayT)
		y <- h*10/16; plot1Row(y,p,11,colA,colB,colP,displayT)
		y <- h*9/16; plot1Row(y,p,13,colA,colB,colP,displayT)
		y <- h*8/16; plot1Row(y,p,15,colA,colB,colP,displayT)

		y <- h*7/16; plot1Row(y,p,17,colA,colB,colP,displayT)
		y <- h*6/16; plot1Row(y,p,19,colA,colB,colP,displayT)
		y <- h*5/16; plot1Row(y,p,21,colA,colB,colP,displayT)
		y <- h*4/16; plot1Row(y,p,23,colA,colB,colP,displayT)
		
		y <- h*3/16; plot1Row(y,p,25,colA,colB,colP,displayT)
		y <- h*2/16; plot1Row(y,p,27,colA,colB,colP,displayT)
		y <- h*1/16; plot1Row(y,p,29,colA,colB,colP,displayT)
		y <- h*0/16; plot1Row(y,p,31,colA,colB,colP,displayT)
	}
}

plotA_withLines <- function(x, y, p, index, colA, colP, lineDown, displayT)
{
	if ((displayT < p$tA1[index]) || (p$tA2[index] < displayT))
	{
		return ()
	}
	if (p$xA1[index] != 0) {
		if (p$StopA[index]) {
			if (p$A_value[index] != '') {
				text(x,y,p$A_value[index],cex=0.8,col=colP)
			  lines(c(x,x-1.5),c(y-0.2,y-lineDown),col=colP)
			}
		} else {
			text(x,y,p$A_value[index],cex=0.8,col=colA)
			lines(c(x,x-1.5),c(y-0.2,y-lineDown),col=colA)
			lines(c(x,x+1.5),c(y-0.2,y-lineDown),col=colA)
		}
	}
}

plotB_withLines <- function(x, y, p, index, colB, colP, lineDown, displayT)
{
	if ((displayT < p$tB1[index]) || (p$tB2[index] < displayT))
	{
		return ()
	}
	if (p$xB1[index] != 0) {
		if (p$StopB[index]) {
			if (p$B_value[index] != '') {
				text(x,y,p$B_value[index],cex=0.8,col=colP)
				lines(c(x,x+1.5),c(y-0.2,y-lineDown),col=colP)
			}
		} else {
			text(x,y,p$B_value[index],cex=0.8,col=colB)
			lines(c(x,x-1.5),c(y-0.2,y-lineDown),col=colB)
			lines(c(x,x+1.5),c(y-0.2,y-lineDown),col=colB)
		}
	}
}

plot2Rows <- function(y1, y2, p, index, colA, colB, colP, lineDown, displayT)
{
	plotA_withLines(2.5,y1,p,index,colA,colP,lineDown,displayT)
	plotB_withLines(8.5,y1,p,index,colB,colP,lineDown,displayT)
	plotA_withLines(2.5,y2,p,index+1,colA,colP,lineDown,displayT)
	plotB_withLines(8.5,y2,p,index+1,colB,colP,lineDown,displayT)
}

display_prev_values <- function(p, ite, colA, colB, colP, nrow, displayT)
{
	h <- nrow+1   # height
	
	if (ite==1)
	{
		y <- h*0.5+1
		plotA_withLines(2.5,y,p,1,colA,colP,0.8,displayT)
		plotB_withLines(8.5,y,p,1,colB,colP,0.8,displayT)
	} else if (ite==2)
	{
		y1 <- h*0.75+1; y2<- h*0.25+1; plot2Rows(y1,y2,p,1,colA,colB,colP,0.8,displayT)
	} else if (ite==3)
	{
		y1 <- h*0.8+1; y2 <- h*0.6+1; plot2Rows(y1,y2,p,1,colA,colB,colP,0.8,displayT)
		y1 <- h*0.4+1; y2 <- h*0.2+1; plot2Rows(y1,y2,p,3,colA,colB,colP,0.8,displayT)
	} else if (ite==4)
	{
		y1 <- h*0.8+0.6; y2 <- h*0.7+0.6; plot2Rows(y1,y2,p,1,colA,colB,colP,0.4,displayT)
		y1 <- h*0.6+0.6; y2 <- h*0.5+0.6; plot2Rows(y1,y2,p,3,colA,colB,colP,0.4,displayT)
		y1 <- h*0.4+0.6; y2 <- h*0.3+0.6; plot2Rows(y1,y2,p,5,colA,colB,colP,0.4,displayT)
		y1 <- h*0.2+0.6; y2 <- h*0.1+0.6; plot2Rows(y1,y2,p,7,colA,colB,colP,0.4,displayT)
	}
}

getARectangle <- function(myData, ite, child)
{
	index <- 0
	if (ite>1)
	{
		for (i in 1:(ite-1))
		{
			index <- index + 2^(i-1)
		}
	}
	index <- index + child
	
	row1 = myData$xA1[index]
	row2 = myData$xA2[index]
	col1 = myData$yA1[index]
	col2 = myData$yA2[index]
	tim1 = myData$tA1[index]
	tim2 = myData$tA2[index]
	
	Stop = myData$StopA[index]
	
	rect <- list(nX=c(row1,row2),nY=c(col1,col2),nT=c(tim1,tim2),Stop=Stop)
	
	return (rect)
}

getBRectangle <- function(myData, ite, child)
{
	index = 0
	if (ite>1)
	{
		for (i in 1:(ite-1))
		{
			index <- index + 2^(i-1)
		}
	}
	index <- index + child
	
	row1 = myData$xB1[index]
	row2 = myData$xB2[index]
	col1 = myData$yB1[index]
	col2 = myData$yB2[index]
	tim1 = myData$tB1[index]
	tim2 = myData$tB2[index]
	
	Stop = myData$StopB[index]
	
	rect <- list(nX=c(row1,row2),nY=c(col1,col2),nT=c(tim1,tim2),Stop=Stop)
	
	return (rect)
}

plotRectangle <- function(rect, color, nrow, displayT)
{
	m1 <- 0.04  # m1 = horizontal margin
	m2 <- 0.03  # m2 = vertical margin
	
	row1 <- rect$nX[1]; row2 <- rect$nX[2]
	col1 <- rect$nY[1]; col2 <- rect$nY[2]
	tim1 <- rect$nT[1]; tim2 <- rect$nT[2]
	
	if ((displayT < tim1) || (tim2 < displayT))
	{
	  return ()
	}
	
	if (rect$Stop) color <- "green"
	
	if (row1 != 0) 
	{
		# top line:
		x <- c(col1-0.5+m2, col2+0.5-m2); y <- c(nrow-row1+1+0.5-m1,nrow-row1+1+0.5-m1)
		lines(x,y,col=color)
		
		# bottom line:
		x <- c(col1-0.5+m2, col2+0.5-m2); y <- c(nrow-row2+1-0.5+m1,nrow-row2+1-0.5+m1)
		lines(x,y,col=color)
		
		# left line:
		x <- c(col1-0.5+m2, col1-0.5+m2); y <- c(nrow-row1+1+0.5-m1,nrow-row2+1-0.5+m1)
		lines(x,y,col=color)
		
		# right line:
		x <- c(col2+0.5-m2, col2+0.5-m2); y <- c(nrow-row1+1+0.5-m1,nrow-row2+1-0.5+m1)
		lines(x,y,col=color)
	}
}

dispPartitionAndTree <- function(ite, nrow, ncol, alpha, lambda, toyexample, myData, 
																 colA, colB, colP, displayT)
{
	myLayout <- matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE)
	nf <- layout(mat=myLayout, widths=c(4,3), heights=c(1,1), respect=FALSE)
	
	# plot numbers:
	par(mar=c(0,0,5,0))
	plot(c(0,ncol+1), c(0,nrow+1), type="n", xlab="", ylab="", axes=FALSE)
	
	iteStr <- ifelse(ite==1," iteration"," iterations")
	title1 <- paste("Partitioning after ", as.character(ite), iteStr, 
									"\nwith z scores (t=", as.character(displayT), ",", 
									"\nalpha=", as.character(alpha), ", ", 
									"lambda=", as.character(lambda), ")", sep="")
	title(main=title1, font.main=1, family="sans", cex.main=1.6)
	
	for (row in 1:nrow)
		for (col in 1:ncol)
		{
			ch <- round(toyexample[row,col,displayT])
			text(col,nrow-row+1,as.character(ch),cex=1)
		}
	
	# plot boundaries:
	maxColourNum = 2^ite
	colour <- rainbow(maxColourNum)
	
	for (child in 1:2^(ite-1))
	{
		rectA = getARectangle(myData, ite, child)
		# plotRectangle(rectA, colour[child])
		plotRectangle(rectA, colA, nrow, displayT)
		
		rectB = getBRectangle(myData, ite, child)
		# plotRectangle(rectB, colour[maxColourNum-child+1])
		plotRectangle(rectB, colB, nrow, displayT)
	}
	
	par(mar=c(0,0,0,0))
	plot(c(0,ncol+1), c(0,nrow+1), type="n", xlab="", ylab="", axes=FALSE)
	
	if (ite>1)
	{
		# display p-values of previous iteration
		values = getValues(myData, ite-1)
		display_prev_values(values, ite-1, colA, colB, colP, nrow, displayT)
	}
	
	# display p-values of current iteration
	values = getValues(myData, ite)
	display_values(values, ite, colA, colB, colP, nrow, displayT)
}

plot1HRow <- function(x, y, m, p, index, colA, colB, colP, displayT)
{
	plotA2(x,y+m,p,index,colA,colP,displayT)
	plotB2(x,y-m,p,index,colB,colP,displayT)
}

plot1HRowWithValues <- function(x, y, m, p, index, colA, colB, colP, displayT)
{
	plotAWithValues(x,y+m,p,index,colA,colP,displayT)
	plotBWithValues(x,y-m,p,index,colB,colP,displayT)
}

plotALine <- function(x1, y1, x2, y2, p, index, colA, colP, displayT)
{
	if (displayT>0)
	{
		if ((displayT < p$tA1[index]) || (p$tA2[index] < displayT)) return ()
	}
	if (p$xA1[index] != 0) {
		if (p$StopA[index])
			lines(c(x1,x2),c(y1,y2),col=colP)
		else
			lines(c(x1,x2),c(y1,y2),col=colA)
	}
}

plotBLine <- function(x1, y1, x2, y2, p, index, colB, colP, displayT)
{
	if (displayT>0)
	{
		if ((displayT < p$tB1[index]) || (p$tB2[index] < displayT)) return ()
	}
	if (p$xB1[index] != 0) {
		if (p$StopB[index])
			lines(c(x1,x2),c(y1,y2),col=colP)
		else
			lines(c(x1,x2),c(y1,y2),col=colB)
	}
}

displayTreeValues <- function(p, ite, colA, colB, colP, nrow, ncol, displayT, maxIte)
{
	h <- nrow+1   # height
	w <- ncol+1   # width
	
	# plot z-scores / p-values:
	x <- w*ite*0.16-w*0.18
	m <- h/2^(ite+1)
	for (i in 1:2^(ite-1))
	{
		y <- ifelse(i==1,2*m,2*m+4*m*(i-1))
		if (ite != maxIte)
			plot1HRow(x, h-y, m, p, i, colA, colB, colP, displayT)
		else
			plot1HRowWithValues(x, h-y, m, p, i, colA, colB, colP, displayT)
	}
	
	# plot joining lines for the tree:
	if (ite>1)
	{
		ox1 <- w*0.08   # offsets for line
		for (i in 1:2^(ite-1))
		{
			y <- ifelse(i==1,2*m,2*m+4*m*(i-1))
			x1 <- x-w*0.15; x2 <- x
			y2 <- h-y+m; y1 <- y2-m
			plotALine(x1+ox1,y1,x2,y2,p,i,colA,colP,displayT)
			y3 <- h-y-m
			plotBLine(x1+ox1,y1,x2,y3,p,i,colB,colP,displayT)
		}
	}
}

dispWholeTree <- function(ite, nrow, ncol, myData, colA, colB, colP, displayT, ignoreTFlag)
{
	# configuration:
	maxIte <- ifelse(ite>5,5,ite)
	
	myLayout <- matrix(c(1,1), nrow=1, ncol=1, byrow=TRUE)
	nf <- layout(mat=myLayout, widths=c(4,3), heights=c(1,1), respect=FALSE)
	
	par(mar=c(0,0,1,0))
	plot(c(0,ncol+1), c(0,nrow+1), type="n", xlab="", ylab="", axes=FALSE)

	if (ignoreTFlag) {
		displayT <- -1
		title1 <- paste("z-scores of partitions for ", as.character(maxIte), 
										" iterations for all times t", sep="")
	} else
		title1 <- paste("z-scores of partitions for ", as.character(maxIte), " iterations at time t=", 
										as.character(displayT), sep="")
	title(main=title1, font.main=1, family="sans", cex.main=1.1)
	
	for (i in 1:maxIte)
	{
		values = getValues(myData, i)
		displayTreeValues(values, i, colA, colB, colP, nrow, ncol, displayT, maxIte)
	}
}

displayValuesInRows <- function(p, ite, colA, colB, colP, nrow, displayT)
{
	h <- nrow+1   # height
	w <- ncol+1   # width
	
	# plot z-scores (rows:cols) [t=times]:
	x <- w*0.01
	m <- h/35
	y <- 2*m
	for (i in 1:2^(ite-1))
	{
		flag <- plotAWithValues(x,h-y,p,i,colA,colP,displayT)
		if (flag) y <- y+2*m
		flag <- plotBWithValues(x,h-y,p,i,colB,colP,displayT)
		if (flag) y <- y+2*m
	}
}

dispPartitionAndValues <- function(ite, nrow, ncol, alpha, lambda, toyexample, myData, 
											  					 colA, colB, colP, displayT)
{
	myLayout <- matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE)
	nf <- layout(mat=myLayout, widths=c(4,3), heights=c(1,1), respect=FALSE)
	
	# plot numbers:
	par(mar=c(0,0,5,0))
	plot(c(0,ncol+1), c(0,nrow+1), type="n", xlab="", ylab="", axes=FALSE)
	
	iteStr <- ifelse(ite==1," iteration"," iterations")
	title1 <- paste("Partitioning after ", as.character(ite), iteStr, 
									"\nwith z scores (t=", as.character(displayT), ",", 
									"\nalpha=", as.character(alpha), ", ", 
									"lambda=", as.character(lambda), ")", sep="")
	title(main=title1, font.main=1, family="sans", cex.main=1.6)
	
	for (row in 1:nrow)
		for (col in 1:ncol)
		{
			ch <- round(toyexample[row,col,displayT])
			text(col,nrow-row+1,as.character(ch),cex=1)
		}
	
	# plot boundaries:
	maxColourNum = 2^ite
	colour <- rainbow(maxColourNum)
	
	for (child in 1:2^(ite-1))
	{
		rectA = getARectangle(myData, ite, child)
		# plotRectangle(rectA, colour[child])
		plotRectangle(rectA, colA, nrow, displayT)
		
		rectB = getBRectangle(myData, ite, child)
		# plotRectangle(rectB, colour[maxColourNum-child+1])
		plotRectangle(rectB, colB, nrow, displayT)
	}
	
	par(mar=c(0,0,5,0))
	plot(c(0,ncol+1), c(0,nrow+1), type="n", xlab="", ylab="", axes=FALSE)

	title1 <- paste("z-score table",
									"\n(row1:row2,col1:col2) ", 
									"\n[t=tim1:tim2]", sep="")
	title(main=title1, font.main=1, family="sans", cex.main=1.6)
	
	# display z-scores of current iteration
	values = getValues(myData, ite)
	displayValuesInRows(values, ite, colA, colB, colP, nrow, displayT)
}

generate3DPoints <- function(lambda, outbreak, xmin, xmax, ymin, ymax, tmin, tmax) {
	area = (xmax - xmin)*(ymax - ymin)
	totalMean = area*lambda
	
	coord.background = NULL
	coord.outbreak = NULL
	for (t in tmin:tmax)
	{
		totalPoints = rpois(1, totalMean)
		
		temp = cbind(runif(totalPoints,xmin,xmax), runif(totalPoints,ymin,ymax), rep(t,totalPoints))
		coord.background = rbind(coord.background, temp)
		
		if(outbreak$nout > 0) {
			for (i in 2:(outbreak$nout+1))
			{
				if ((outbreak[[i]]$tmin <= t) && (t <= outbreak[[i]]$tmax)) {
					xvector = runif(outbreak[[i]]$ncounts,outbreak[[i]]$xmin,outbreak[[i]]$xmax)
					yvector = runif(outbreak[[i]]$ncounts,outbreak[[i]]$ymin,outbreak[[i]]$ymax)
					tvector = rep(t,outbreak[[i]]$ncounts)
					coord.xyt = cbind(xvector, yvector, tvector)
					coord.outbreak = rbind(coord.outbreak, coord.xyt)
				}
			}
		}
	}
	
	return (list(obs.all = rbind(coord.background,coord.outbreak), 
							 obs.background = coord.background, 
							 obs.outbreak = coord.outbreak))
}

# Generate ncounts outbreak points (x, y, t) for each time slice, 
# where xmin <= x <= xmax (row), ymin <= y <= ymax (column), and tmin <= t <= tmax (time):
#
# total number of points generated = ncounts * (tmax-tmin+1)
#
generateOutbreakList <- function(ncounts, xmin, xmax, ymin, ymax, tmin, tmax) {
	return (list(ncounts=ncounts, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, tmin=tmin, tmax=tmax))
}

convertPointsToGrid <- function(points, xmin, xmax, ymin, ymax, tmin, tmax) {
	X = xmax - xmin; Y = ymax - ymin; T = tmax - tmin + 1
	dim = c(X,Y,T)   # dimension of the space (X x Y x T)
	Obs = array(rep(0,prod(dim)), dim) # 3D observation data
	
	npoints = dim(points)[1]
	for (i in 1:npoints)
	{
		x = ceiling(points[i,1])
		y = ceiling(points[i,2])
		t = points[i,3]
		Obs[x,y,t] = Obs[x,y,t]+1
	}
	
	return (Obs)
}

plotBoundaries <- function(xmin, xmax, ymin, ymax, colBoundary) {
	plot.new()   # clear plot
	
	myLayout <- matrix(c(1,1), nrow=1, ncol=1, byrow=TRUE)
	nf <- layout(mat=myLayout, widths=c(4,3), heights=c(1,1), respect=FALSE)
	
	par(mar=c(0,0,1,0))
	plot(c(ymin,ymax), c(xmin,xmax), type="n", xlab="", ylab="", axes=FALSE)
	
	# plot boundaries: horizontal lines:
	for (row in xmin:xmax)
	{
		x <- c(ymin, ymax); y <- c(row,row)
		lines(x,y,col=colBoundary)
	}
	
	# plot boundaries: vertical lines:
	for (col in ymin:ymax)
	{
		x <- c(col, col); y <- c(xmin,xmax)
		lines(x,y,col=colBoundary)
	}
}

plot3DPoints <- function(points, Obs, xmin, xmax, ymin, ymax, tmin, tmax, colBoundary, colText) {
	npoints = dim(points)[1]
	
	for (t in tmin:tmax)
	{
		plotBoundaries(xmin, xmax, ymin, ymax, colBoundary)
		
		title1 <- paste("Points & counts for t=", as.character(t), sep="")
		title(main=title1, font.main=1, family="sans", cex.main=1.6)
		
		# plot text
		Obs2 <- Obs[,,t]
		for (row in 1:xmax)
			for (col in 1:ymax)
			{
				ch <- Obs2[row,col]
				text(col-0.5,xmax-row+0.5,as.character(ch),cex=1,col=colText)
			}
		
		# plot points
		current_x <- NULL; current_y <- NULL
		for (i in 1:npoints)
		{
			current_t = points[i,3]
			if (t == current_t)
			{
				current_x <- c(current_x, points[i,1])
				current_y <- c(current_y, points[i,2])
			}
		}
		#points(current_y, xmax-current_x, pch='.', cex=2)
		points(current_y, xmax-current_x, pch=20, cex=0.75)
		
		ch = readline(prompt="Press [enter] to continue, or [t] to terminate plotting of points: ")
		if ((ch == "t") || (ch == 'T')) break
	}
}

readAegissData <- function(filename, T1, T2) {
	data <- read.table(filename, header=T)
	
	xyAll <- cbind(data$x, data$y, data$t)
	
	xmin = min(data$x); xmax = max(data$x)
	ymin = min(data$y); ymax = max(data$y)
	
	m <- dim(data)[1]
	
	vecx1 <- rep(0,m); vecy1 <- rep(0,m)    # before
	vecx2 <- rep(0,m); vecy2 <- rep(0,m)    # after
	for (i in 1:m)
	{
		if (T1 <= data$t[i]) {
			if (data$t[i] <= T2) {
				vecx1[i] <- data$x[i]; vecy1[i] <- data$y[i]
			}
			else {
				vecx2[i] <- data$x[i]; vecy2[i] <- data$y[i]
			}
		}
	}
	xyBefore <- cbind(vecx1, vecy1)
	xyAfter <- cbind(vecx2, vecy2)
	
	# override (min, max) on x- and y-axes:
	xmin = 400000; xmax = 490000
	ymin = 90000;  ymax = 170000
	xyDim = list(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
	
	return (list(All=xyAll, Before=xyBefore, After=xyAfter, Dim=xyDim))
}

readAegissPoly <- function(filename) {
	vec1 <- scan(filename, what=character())
	n <- length(vec1)
	vec2 <- rep(0,n)
	for (i in 1:n)
		vec2[i] <- as.numeric(vec1[i])
	
	poly <- matrix(vec2, 120, 2, T)
	return (poly)
}

plotAegiss <- function(xyPart, poly, title, pauseFlag) {
	par(pty="s")      # generate a square plotting region
	#par(pty="m")      # generate a maximal plotting region
	polymap(poly)
	pointmap(xyPart, add=T, pch=".")
	#pointmap(xyPart, add=T, pch=19, cex=0.25)
	
	title(main=title, font.main=1, family="sans", cex.main=1.6)
	
	if (pauseFlag)
		invisible(readline(prompt="Press [enter] to continue"))
}

plotAegissBoundaries <- function(poly) {
	n = dim(poly)[1]
	x0 = x1 = poly[1,1]
	y0 = y1 = poly[1,2]
	lineWidth = 1   # 1 or 2
	for (i in 2:n) {
		x2 = poly[i,1]; y2 = poly[i,2]
		lines(c(x1,x2),c(y1,y2),lwd=lineWidth)
		x1 = x2; y1 = y2;
	}
	lines(c(x1,x0),c(y1,y0),lwd=lineWidth)
}

plotAegissPoints <- function(xyPart, colP) {
# 	npoints = dim(xyPart)[1]
# 	current_x <- NULL; current_y <- NULL
# 	for (i in 1:npoints)
# 	{
# 		current_x <- c(current_x, xyPart[i,1])
# 		current_y <- c(current_y, xyPart[i,2])
# 	}
	#points(current_x, current_y, pch='.', cex=2)
	points(xyPart[,1], xyPart[,2], pch='.', cex=0.75, col=colP)
	#points(current_x, current_y, pch=20, cex=0.75)
}

plotAegissGrid <- function(Dim, Obs, nY, nX)
{
	xmin = Dim$xmin; xmax = Dim$xmax; xlen = (xmax-xmin)/nX
	ymin = Dim$ymin; ymax = Dim$ymax; ylen = (ymax-ymin)/nY

	# plot horizontal lines:
	x1 = xmin; x2 = xmax
	for (i in 0:nY) {
		y1 = y2 = ymin+i*ylen
		lines(c(x1,x2),c(y1,y2),col='red')
	}

	# plot vertical lines:
	y1 = ymin; y2 = ymax
	for (j in 0:nX) {
		x1 = x2 = xmin+j*xlen
		lines(c(x1,x2),c(y1,y2),col='red')
	}
	
	# plot text for counts:
	#text(x,y,p$A_value[index],cex=0.8,col=colorA)
	for (row in 1:nY)
		for (col in 1:nX)
		{
			ch <- round(Obs[row,col])
			if (ch>0) {
				y = ymin + ylen*(nY-row+1-0.5)
				x = xmin + xlen*(col-0.5)
				if (nX < 15)
					text(x,y,as.character(ch),cex=1,col='black')
				else
					text(x,y,as.character(ch),cex=0.7,col='black')
			}
		}
}

convertAegissPointsToGrid <- function(xy, nY, nX, T2, nD) {
	xmin = xy$Dim$xmin; xmax = xy$Dim$xmax; xlen = (xmax-xmin)/nX
	ymin = xy$Dim$ymin; ymax = xy$Dim$ymax; ylen = (ymax-ymin)/nY
	
	dim = c(nY,nX)   # dimension of the space (nY x nX)
	ObsAll = array(rep(0,prod(dim)), dim) # 2D observation data
	ObsBefore = array(rep(0,prod(dim)), dim) # 2D observation data
	ObsAfter = array(rep(0,prod(dim)), dim) # 2D observation data
	
	npoints = dim(xy$All)[1]
	for (i in 1:npoints)
	{
		x = xy$All[i,1]
		y = xy$All[i,2]
		ix = ceiling((x-xmin)/xlen)
		iy = nY-ceiling((y-ymin)/ylen)+1
		ObsAll[iy,ix] = ObsAll[iy,ix]+1
	}

	for (i in 1:npoints)
	{
		x = xy$Before[i,1]
		y = xy$Before[i,2]
		if ((x <= 0) || (y <= 0)) {
			break
		}
		ix = ceiling((x-xmin)/xlen)
		iy = nY-ceiling((y-ymin)/ylen)+1
		ObsBefore[iy,ix] = ObsBefore[iy,ix]+1
	}

	for (i in 1:npoints)
	{
		x = xy$After[i,1]
		y = xy$After[i,2]
		if ((x <= 0) || (y <= 0)) {
			next
		}
		ix = ceiling((x-xmin)/xlen)
		iy = nY-ceiling((y-ymin)/ylen)+1
		ObsAfter[iy,ix] = ObsAfter[iy,ix]+1
	}

	nT <- ceiling(max(xy$All[,3])/nD)   # number of time slices (nD days per slice)
	dim2 = c(nY,nX,nT)                  # dimension of the space (nY x nX)
	ObsCounts = array(rep(0,prod(dim2)), dim2) # 3D observation data
	for (i in 1:npoints)
	{
		x = xy$All[i,1]
		y = xy$All[i,2]
		t = xy$All[i,3]
		ix = ceiling((x-xmin)/xlen)
		iy = nY-ceiling((y-ymin)/ylen)+1
		it = ceiling(t/nD)
		ObsCounts[iy,ix,it] = ObsCounts[iy,ix,it]+1
	}
	
	# checking counts per nD days:
	nDays=1
	maxD = ceiling(max(xy$All[,3])/nDays)
	Obs1D = rep(0,maxD)
	for (i in 1:npoints)
	{
		t = xy$All[i,3]
		it = ceiling(t/nDays)
		Obs1D[it] = Obs1D[it]+1
	}
	sum = sum(Obs1D)

	nDays=7
	maxD = ceiling(max(xy$All[,3])/nDays)
	Obs7D = rep(0,maxD)
	for (i in 1:npoints)
	{
		t = xy$All[i,3]
		it = ceiling(t/nDays)
		Obs7D[it] = Obs7D[it]+1
	}

	nDays=30
	maxD = ceiling(max(xy$All[,3])/nDays)
	Obs30D = rep(0,maxD)
	for (i in 1:npoints)
	{
		t = xy$All[i,3]
		it = ceiling(t/nDays)
		Obs30D[it] = Obs30D[it]+1
	}
	
# 	for (i in 1:nY)
# 		for (j in 1:nX)
# 			ObsMeans[i,j,t] = ObsBefore[i,j] / T2 * nT
	
	# checking sums:
# 	Obs2 = array(rep(0,prod(dim)), dim) # 2D observation data
# 	for (i in 1:nY)
# 		for (j in 1:nX)
# 			Obs2[i,j] = ObsBefore[i,j] + ObsAfter[i,j] - ObsAll[i,j]
# 	sumObs = sum(Obs2)
	
	return (list(All=ObsAll, Before=ObsBefore, After=ObsAfter, Counts=ObsCounts, 
							 Days1=Obs1D, Days7=Obs7D, Days30=Obs30D))
}

plotAegiss2 <- function(xyPart, poly, Dim, title, Obs, nY, nX, colP, pauseFlag) {
	xmin = Dim$xmin; xmax = Dim$xmax
	ymin = Dim$ymin;  ymax = Dim$ymax
	
	myLayout <- matrix(c(1,1), nrow=1, ncol=1, byrow=TRUE)
	nf <- layout(mat=myLayout)
	
	par(mar=c(4,2.5,2,2))
	plot(c(xmin,xmax), c(ymin,ymax), type="n", xlab="", ylab="", axes=TRUE)

	plotAegissBoundaries(poly)
	plotAegissPoints(xyPart,colP)
	
	title(main=title, font.main=1, family="sans", cex.main=1.6)
	
	plotAegissGrid(Dim, Obs$All, nY, nX)
	
	if (pauseFlag)
		invisible(readline(prompt="Press [enter] to continue"))
}

plotDaySeries2 <- function(timeSeries, nDays) {
	nD2 <- length(timeSeries)
	x = seq(1,nD2)
	
	plot(timeSeries,pch=3,ylim=range(0,max(timeSeries)))
	lines(x,timeSeries,lty=1)
	
	# draw 1st year line:
	xT2 = ceiling(365/nDays)
	lines(c(xT2,xT2),c(0,max(timeSeries)),lty=2,col='red')
	
	# draw 2nd year line:
	xT2 = ceiling(365*2/nDays)
	lines(c(xT2,xT2),c(0,max(timeSeries)),lty=2,col='red')
}

plotDaySeries <- function(Obs) {
	myLayout <- matrix(c(1,2,3), nrow=3, ncol=1, byrow=TRUE)
	nf <- layout(mat=myLayout, widths=c(4,3), heights=c(1,1), respect=FALSE)
	
	plotDaySeries2(Obs$Days1, 1)
	title1 <- "Time series for all regions for each day"
	title(main=title1, font.main=1, family="sans", cex.main=1.6)

	plotDaySeries2(Obs$Days7, 7)
	title1 <- "For each week"
	title(main=title1, font.main=1, family="sans", cex.main=1.6)
	
	plotDaySeries2(Obs$Days30, 30)
	title1 <- "For each 30 days"
	title(main=title1, font.main=1, family="sans", cex.main=1.6)
}

plotTimeSeries <- function(timeSeries, T2, nD, cellIndex, means) {
	nT <- length(timeSeries)
	x = seq(1,nT)
	ymean = means
	
	plot(timeSeries,pch=3,ylim=range(0,max(timeSeries)))
	lines(x,timeSeries,lty=1)

	# means for 1st year:
	points(x[1:13],ymean[1:13],pch=3,col='red')
	lines(x[1:13],ymean[1:13],lty=2,col='red')

	# means for 2nd year:
	points(x[13:25],ymean[13:25],pch=3,col='blue')
	lines(x[13:25],ymean[13:25],lty=2,col='blue')
	
	# means for 3rd year:
	points(x[25:nT],ymean[25:nT],pch=3,col='red')
	lines(x[25:nT],ymean[25:nT],lty=2,col='red')
	
	# draw 1st year line:
	xT2 = ceiling(365/nD)
	lines(c(xT2,xT2),c(0,max(timeSeries)),lty=2,col='green')
	
	# draw 2nd year line:
	#xT2 = ceiling(T2/nD)
	xT2 = ceiling(365*2/nD)
	lines(c(xT2,xT2),c(0,max(timeSeries)),lty=2,col='green')
	
	row = cellIndex[1]
	col = cellIndex[2]
	title1 <- paste("Time series for (", as.character(row), ",", as.character(col), ") with sum = ",
									as.character(sum(timeSeries)), sep="")
	title(main=title1, font.main=1, family="sans", cex.main=1.6)
	
	legend("topright",c("counts","mean"),lty=c(1,2),col=c("black","red"),text.col=c("black","red"))
	
	str = paste("Time series for (", as.character(row), ",", as.character(col), ")", sep="")
	print(str)
	print(timeSeries)
	
	str = paste("total sum = ",as.character(sum(timeSeries)), sep="")
	print(str)
	
	it = ceiling(T2/nD)
	mean = 0
	for (t in 1:it)
		mean = mean + timeSeries[t]
	mean = mean/it
	str = paste("mean = ",as.character(mean), sep="")
	print(str)
}

findTimeSeries <- function(Obs, cellIndices, nT, nD) {
	numIndices <- dim(cellIndices)[1]
	
	# find sums & means:
	it = ceiling(T2/nD)
	
	sums <- rep(0,numIndices)
	
	dimT = c(numIndices,nT)
	means <- array(rep(0,prod(dimT)), dimT)
	timeSeries <- array(rep(0,prod(dimT)), dimT)
	
	for (i in 1:numIndices) {
		row = cellIndices[i,1]
		col = cellIndices[i,2]
		
		sum = 0
		for (t in 1:nT)
			sum = sum + Obs$Counts[row,col,t]
		sums[i] = sum

		# find average of each month from 2 year data:
# 		for (t in 1:12) {
# 			means[i,t] = (Obs$Counts[row,col,t]+Obs$Counts[row,col,t+12])/2
# 		}
# 		for (t in 13:nT) {
# 			means[i,t] = means[i,t-12]
# 		}
    maxt = floor(365/nD)
		for (t in 1:maxt) {
			means[i,t] = (Obs$Counts[row,col,t]+Obs$Counts[row,col,t+maxt])/2
		}
		for (t in (maxt+1):nT) {
			means[i,t] = means[i,t-maxt]
		}
		
		timeSeries[i,] <- Obs$Counts[row,col,]
	}
	
	# row index, col index, sum, mean, nT time series
	return (list(indices=cellIndices, sums=sums, means=means, timeSeries=timeSeries))
}

findOutbreakData <- function(Obs, T2, nD) {
	nY <- dim(Obs$Counts)[1]
	nX <- dim(Obs$Counts)[2]
	nT <- dim(Obs$Counts)[3]
	
	it = ceiling(T2/nD)
	it2 = nT - it
	
	# find data and means for recursive partitioning:
	dimD = c(nY,nX,it2)
	ObsData <- array(rep(0,prod(dimD)), dimD)
	ObsMeans <- array(rep(0,prod(dimD)), dimD)
	
	num1Period = 12
	
	for (y in 1:nY)
		for (x in 1:nX) {
			ObsData[y,x,] = Obs$Counts[y,x,(it+1):nT]
			
			for (i in (it+1):nT) {
				k <- i %% num1Period
				if (k == 0) k = num1Period
				mean = (Obs$Counts[y,x,k] + Obs$Counts[y,x,k+num1Period])/2
				ObsMeans[y,x,i-it] = mean
			}
		}
	
	# find data and means at time it for temporal smoothing:
	dimD = c(nY,nX)
	ObsData0 <- array(rep(0,prod(dimD)), dimD)
	ObsMeans0 <- array(rep(0,prod(dimD)), dimD)
	
	for (y in 1:nY)
		for (x in 1:nX) {
			ObsData0[y,x] = Obs$Counts[y,x,it]
			ObsMeans0[y,x] = Obs$Counts[y,x,it]
		}
	
	return (list(data=ObsData, means=ObsMeans, data0=ObsData0, means0=ObsMeans0))
}

plotTimeSlice <- function(xmin, xmax, ymin, ymax, xyAll, poly, nD, it, current_t, colP, 
													Dim, data, nY, nX) {
	par(mar=c(4,2.5,2,2))
	plot(c(xmin,xmax), c(ymin,ymax), type="n", xlab="", ylab="", axes=TRUE)
	
	plotAegissBoundaries(poly)
	
	npoints = dim(xyAll)[1]
	for (i in 1:npoints)
	{
		t = ceiling(xy$All[i,3]/nD)
		if (t == (it+current_t)) {
			x = xy$All[i,1]
			y = xy$All[i,2]
			points(x, y, pch=16, cex=0.75, col=colP)
		}
	}
	
	plotAegissGrid(Dim, data, nY, nX)
	
}

plotOutbreakTimeSlices <- function(xyAll, poly, Dim, colP, T2, nD, nY, nX, 
																	 data, means, original, nrow, ncol, current_t) {
	myLayout <- matrix(c(1,2,3,4), nrow=2, ncol=2, byrow=TRUE)
	nf <- layout(mat=myLayout, widths=c(1,1), heights=c(1,1), respect=FALSE)
	layout.show(n = 4)
	
	xmin = Dim$xmin; xmax = Dim$xmax
	ymin = Dim$ymin;  ymax = Dim$ymax

	it = ceiling(T2/nD)

	# plot 1st graph - data:
	plotTimeSlice(xmin, xmax, ymin, ymax, xyAll, poly, nD, it, current_t, colP, 
								Dim, data, nY, nX)

	title1 <- paste("Data for time slice t=", as.character(current_t), sep="")
	title(main=title1, font.main=1, family="sans", cex.main=1.6)

	# plot 2nd graph - means:
	plotTimeSlice(xmin, xmax, ymin, ymax, xyAll, poly, nD, it, current_t, colP, 
								Dim, means, nY, nX)
	
	title(main="Means", font.main=1, family="sans", cex.main=1.6)

	# plot 3rd graph - data-means>0:
	diff = data - means
	diff[diff<0] = 0
	plotTimeSlice(xmin, xmax, ymin, ymax, xyAll, poly, nD, it, current_t, colP, 
								Dim, diff, nY, nX)
	
	title(main="Data - Means > 0", font.main=1, family="sans", cex.main=1.6)

	# plot 4th graph - original:
	plotTimeSlice(xmin, xmax, ymin, ymax, xyAll, poly, nD, it, current_t, colP, 
								Dim, original, nY, nX)
	
	title(main="Original (before smoothing)", font.main=1, family="sans", cex.main=1.6)
}

