# This file works with Outbreak_Detection2.R & Outbreak_Detection3.R

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

####

# inputs:
#       scounts: the current matrix O2 to be analysed for row/col partitioning
#       smus: corresponding matrix with expected counts (same dimension as scounts)
#       gen.dim: $row: initial row and last row indices
#                $column: initial col and last col indices
#
# choosing the best row/col partition A for 1 iteration, and returning part: 
#       col.part = FALSE if choosing row partition
#                  TRUE if choosing col partition
#       partition = [1] = initial index for row/col partition A
#                   [2] = last index for row/col partition
#       gen.dim: $row / $column updated with the chosen partition A
#
recursive.partition <- function(criteria,scounts,smus,gen.dim) {
	if (criteria==1)
		part <- recursive.partition.p(scounts,smus,gen.dim)
	else
		part <- recursive.partition.z(scounts,smus,gen.dim)
	
	return (part)
}

# partitioning by p-values:
recursive.partition.p <- function(scounts,smus,gen.dim) {
	# yr, yrmu: row sums
	# yc, ycmu: column sums
  if (gen.dim$row[1]!=gen.dim$row[2] & gen.dim$column[1]!=gen.dim$column[2]){
    yr <- apply(scounts,1,sum); yrmu<-apply(smus,1,sum)
    yc <- apply(scounts,2,sum); ycmu<-apply(smus,2,sum)
  } else {
    if (gen.dim$column[1]!=gen.dim$column[2]){
    	yr <- sum(scounts); yrmu<-sum(smus)
      yc <- scounts; ycmu<-smus
  	} else {
  		yr <- scounts; yrmu<-smus; 
    	yc <- sum(scounts); ycmu<-sum(smus)
  	}
  }
  
  nrow <- length(yr)
  ncol <- length(yc)
  
  # pr.up: p-values in upward row direction, pr.down: p-values in downward row direction
  # pc.up: p-values in rightward col direction, pc.down: p-values in leftward col direction
  pr.up   <- 1-ppois(cumsum(yr), cumsum(yrmu))
  pr.down <- 1-ppois(cumsum(yr[nrow:1]), cumsum(yrmu[nrow:1]))
  pc.up   <- 1-ppois(cumsum(yc), cumsum(ycmu))
  pc.down <- 1-ppois(cumsum(yc[ncol:1]), cumsum(ycmu[ncol:1]))
  
  # make last p-value = 1 in order to eliminate the case of no partitioning:
  pr.up[nrow] = 1; pr.down[nrow] = 1
  pc.up[ncol] = 1; pc.down[ncol] = 1
  
  # tem.row: min p-values for row partition
  # tem.col: min p-values for col partition
  tem.row <- c(min(pr.up),min(pr.down))
  tem.col <- c(min(pc.up),min(pc.down))
  
  # Look for the min p-value for row sums - summing over rows up and then down
  # p1: index for min p-value for row partition in upward direction
  # p2: index for min p-value for row partition in downward direction
  p1 <- c(1:nrow)[pr.up==tem.row[1]]
  p2 <- c(nrow:1)[pr.down==tem.row[2]]
  
  # Look for the min p-value for column sums - summing over columns up and then down
  # p3: index for min p-value for col partition in upward direction
  # p4: index for min p-value for col partition in downward direction
  p3<-c(1:ncol)[pc.up==tem.col[1]]
  p4<-c(ncol:1)[pc.down==tem.col[2]]
  
  # in case there are multiple min p-values
  if (length(p1)>1) p1 <- p1[1]
  if (length(p2)>1) p2 <- p2[1]
  if (length(p3)>1) p3 <- p3[1]
  if (length(p4)>1) p4 <- p4[1]
  
  # Recording all potential "best" partitions
  # potential2: 1st row index, row index for max upward z-score, 
  #                            row index for max downward z-score, last row index,
  #             1st col index, col index for max leftward z-score,
  #                            col index for max rightward z-score, last col index.
  potential2 <- matrix(c(gen.dim$row[1],gen.dim$row[1]+p1-1,
  											 gen.dim$row[1]+p2-1,gen.dim$row[2],
  											 gen.dim$column[1],gen.dim$column[1]+p3-1,
  											 gen.dim$column[1]+p4-1,gen.dim$column[2]),4,2,byrow=T)
  
  # Check whether no partitioning is necessary of the rows
  trmin <- min(tem.row)  # overall min p-value in row partition
  if (gen.dim$row[1]==gen.dim$row[2]) {
  	row1 <- 1
  } else {
    if(tem.row[1]==trmin) {
    	# row1 = 1 to signify best row partition in upward direction
      row1 <- 1
    } else {
    	# row1 = 2 to signify best row partition in downward direction
      row1 <- 2
    }
  } 
  
  # Check whether no partitioning is necessary of the columns
  tcmin <- min(tem.col)  # overall max z-score in col partition
  if(gen.dim$column[1]==gen.dim$column[2]){
    row2 <- 3
  } else {
    if(tem.col[1]==tcmin){
    	# row2 = 3 to signify best col partition in leftward direction
      row2 <- 3
    } else {
    	# row2 = 4 to signify best col partition in rightward direction
      row2 <- 4
    }
  } 
  
  # tmax: overall max normalized z-score for both row & col partitions
  # ind: 1 to signify best row partition in upward direction
  #      2 to signify best row partition in downward direction
  #      3 to signify best col partition in leftward direction
  #      4 to signify best col partition in rightward direction
  # col.part: FALSE if choosing row partition
  #           TRUE if choosing col partition
  tmin = min(tem.row, tem.col)
  if(trmin==tmin) {
  	ind<-row1
    col.part<-FALSE
  	if (ind==1) {
  		A_value = pr.up[p1]
  		B_value = pr.down[nrow-p1]
  	} else {
  		p2a <- c(1:nrow)[pr.down==tmin]
  		if (length(p2a)>1) p2a <- p2a[1]

  		A_value = pr.down[p2a]
  		B_value = pr.up[nrow-p2a]
  	}
  } else {
  	ind<-row2
    col.part<-TRUE
  	if (ind==3) {
  		A_value = pc.up[p3]
  		B_value = pc.down[ncol-p3]
  	} else {
  		p4a<-c(1:ncol)[pc.down==tmin]
  		if (length(p4a)>1) p4a <- p4a[1]

  		A_value = pc.down[p4a]
  		B_value = pc.up[ncol-p4a]
  	}
  }
  
  # partition: [1] = initial index for row/col partition
  #            [2] = last index for row/col partition
  partition<-potential2[ind,]
  
  if (col.part==T) 
  	gen.dim$column<-partition 
  else 
  	gen.dim$row<-partition
  
  values <- list(A_value=A_value, B_value=B_value)
 
  # part: col.part = FALSE if choosing row partition
  #                  TRUE if choosing col partition
  #       partition = [1] = initial index for row/col partition
  #                   [2] = last index for row/col partition
  #       gen.dim: $row & $column updated with the chosen partition
  part <- list(col.part=col.part,partition=partition,gen.dim=gen.dim,values=values)
  
  part
}

# partitioning by z-scores:
recursive.partition.z <- function(scounts,smus,gen.dim) {
	# yr, yrmu: row sums
	# yc, ycmu: column sums
	if (gen.dim$row[1]!=gen.dim$row[2] & gen.dim$column[1]!=gen.dim$column[2]){
		yr <- apply(scounts,1,sum); yrmu<-apply(smus,1,sum)
		yc <- apply(scounts,2,sum); ycmu<-apply(smus,2,sum)
	} else {
		if (gen.dim$column[1]!=gen.dim$column[2]){
			yr <- sum(scounts); yrmu<-sum(smus)
			yc <- scounts; ycmu <- smus
		} else {
			yr <- scounts; yrmu<-smus
			yc <- sum(scounts); ycmu<-sum(smus)
		}
	}
	
	nrow <- length(yr)
	ncol <- length(yc)
	
	# pr.up: z-scores in upward row direction, pr.down: z-scores in downward row direction
	# pc.up: z-scores in rightward col direction, pc.down: z-scores in leftward col direction
	pr.up   <- 2*(sqrt(cumsum(yr))-sqrt(cumsum(yrmu)))
	pr.down <- 2*(sqrt(cumsum(yr[length(yr):1]))-sqrt(cumsum(yrmu[length(yr):1])))
	pc.up   <- 2*(sqrt(cumsum(yc))-sqrt(cumsum(ycmu)))
	pc.down <- 2*(sqrt(cumsum(yc[length(yc):1]))-sqrt(cumsum(ycmu[length(yc):1])))
	
	# make last z-score = min. z-score - 1 in order to eliminate the case of no partitioning:
	zmin = min(c(pr.up, pr.down, pc.up, pc.down))
	pr.up[nrow] = zmin-1; pr.down[nrow] = zmin-1
	pc.up[ncol] = zmin-1; pc.down[ncol] = zmin-1
	
	# tem.row: max z-scores for row partition
	# tem.col: max z-scores for col partition
	tem.row <- c(max(pr.up),max(pr.down))
	tem.col <- c(max(pc.up),max(pc.down))
	
	# Look for the max z-score for row sums - summing over rows up and then down
	# p1: index for max z-score for row partition in upward direction
	# p2: index for max z-score for row partition in downward direction
	p1 <- c(1:nrow)[pr.up==tem.row[1]]
	p2 <- c(nrow:1)[pr.down==tem.row[2]]
	
	# Look for the max z-score for column sums - summing over columns up and then down
	# p3: index for max z-score for col partition in upward direction
	# p4: index for max z-score for col partition in downward direction
	p3<-c(1:ncol)[pc.up==tem.col[1]]
	p4<-c(ncol:1)[pc.down==tem.col[2]]
	
	# sometimes, there are multiple max z-scores
	if (length(p1)>1) p1 <- p1[1]
	if (length(p2)>1) p2 <- p2[1]
	if (length(p3)>1) p3 <- p3[1]
	if (length(p4)>1) p4 <- p4[1]
	
	# Recording all potential "best" partitions
	# potential2: 1st row index, row index for max upward z-score, 
	#                            row index for max downward z-score, last row index,
	#             1st col index, col index for max leftward z-score,
	#                            col index for max rightward z-score, last col index.
	potential2 <- matrix(c(gen.dim$row[1],gen.dim$row[1]+p1-1,
												 gen.dim$row[1]+p2-1,gen.dim$row[2],
												 gen.dim$column[1],gen.dim$column[1]+p3-1,
												 gen.dim$column[1]+p4-1,gen.dim$column[2]),4,2,byrow=T)
	
	# Check whether no partitioning is necessary of the rows
	trmax <- max(tem.row)  # overall max z-score in row partition
	if (gen.dim$row[1]==gen.dim$row[2]) {
		row1 <- 1
	} else {
		if(tem.row[1]==trmax) {
			# row1 = 1 to signify best row partition in upward direction
			row1 <- 1
		} else {
			# row1 = 2 to signify best row partition in downward direction
			row1 <- 2
		}
	} 
	
	# Check whether no partitioning is necessary of the columns
	tcmax <- max(tem.col)  # overall max z-score in col partition
	if(gen.dim$column[1]==gen.dim$column[2]){
		row2 <- 3
	} else {
		if(tem.col[1]==tcmax){
			# row2 = 3 to signify best col partition in leftward direction
			row2 <- 3
		} else {
			# row2 = 4 to signify best col partition in rightward direction
			row2 <- 4
		}
	} 
	
	# tmax: overall max normalized z-score for both row & col partitions
	# ind: 1 to signify best row partition in upward direction
	#      2 to signify best row partition in downward direction
	#      3 to signify best col partition in leftward direction
	#      4 to signify best col partition in rightward direction
	# col.part: FALSE if choosing row partition
	#           TRUE if choosing col partition
	tmax = max(tem.row, tem.col)
	if(trmax==tmax) {
		ind<-row1
		col.part<-FALSE
		if (ind==1) {
			A_value = pr.up[p1]
			B_value = pr.down[nrow-p1]
		} else {
			p2a <- c(1:nrow)[pr.down==tmax]
			if (length(p2a)>1) p2a <- p2a[1]

			A_value = pr.down[p2a]
			B_value = pr.up[nrow-p2a]
		}
	} else {
		ind<-row2
		col.part<-TRUE
		if (ind==3) {
			A_value = pc.up[p3]
			B_value = pc.down[ncol-p3]
		} else {
			p4a<-c(1:ncol)[pc.down==tmax]
			if (length(p4a)>1) p4a <- p4a[1]
			
			A_value = pc.down[p4a]
			B_value = pc.up[ncol-p4a]
		}
	}
	
	# partition: [1] = initial index for row/col partition
	#            [2] = last index for row/col partition
	partition<-potential2[ind,]
	
	if (col.part==T) 
		gen.dim$column<-partition 
	else 
		gen.dim$row<-partition
	
	values <- list(A_value=A_value, B_value=B_value)
	
	# part: col.part = FALSE if choosing row partition
	#                  TRUE if choosing col partition
	#       partition = [1] = initial index for row/col partition
	#                   [2] = last index for row/col partition
	#       gen.dim: $row & $column updated with the chosen partition
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
  
  if(col.part==TRUE) {
    # hp.count = extracted desired col partition
    # hp.mu = similar to hp.count of the same corresponding dimension
    # gen.dim1 = extracted row & col dimensions
    hp.count <- scount[gen.dim$row[1]:gen.dim$row[2],partition[1]:partition[2]];
    hp.mu <- smu[gen.dim$row[1]:gen.dim$row[2],partition[1]:partition[2]]
    gen.dim1$column<-partition
    
    # gen.dim2 = row & col dimensions of the other col partition
    if (partition[1]==gen.dim$column[1]) {
      if (partition[2]!=gen.dim$column[2])
      	gen.dim2$column<-c(partition[2]+1,gen.dim$column[2])
      else
      	# should not arrive here anyway
      	gen.dim2$column<-c(gen.dim$column[2],gen.dim$column[2])
    } else
    	gen.dim2$column<-c(gen.dim$column[1],partition[1]-1)
    
    # lp.count = the other col partition
    # lp.mu = similar to lp.count of the same corresponding dimension
    lp.count <- scount[gen.dim2$row[1]:gen.dim2$row[2],gen.dim2$column[1]:gen.dim2$column[2]];
    lp.mu <- smu[gen.dim2$row[1]:gen.dim2$row[2],gen.dim2$column[1]:gen.dim2$column[2]] 
  } else {
  	# hp.count = extracted desired row partition
  	# hp.mu = similar to hp.count of the same corresponding dimension
  	# gen.dim1 = extracted row & col dimensions
    hp.count <- scount[c(partition[1]:partition[2]),gen.dim$column[1]:gen.dim$column[2]];
    hp.mu <- smu[c(partition[1]:partition[2]),gen.dim$column[1]:gen.dim$column[2]]
    gen.dim1$row <- partition
    
  	# gen.dim2 = row & col dimensions of the other row partition
  	if (partition[1]==gen.dim$row[1]) {
      if (partition[2]!=gen.dim$row[2])
      	gen.dim2$row<-c(partition[2]+1,gen.dim$row[2])
      else
      	# should not arrive here anyway
      	gen.dim2$row<-c(gen.dim$row[2],gen.dim$row[2])
    } else
    	gen.dim2$row<-c(gen.dim$row[1],partition[1]-1)
    
  	# lp.count = the other row partition
  	# lp.mu = similar to lp.count of the same corresponding dimension
  	lp.count<-scount[gen.dim2$row[1]:gen.dim2$row[2],gen.dim$column[1]:gen.dim$column[2]];
    lp.mu<-smu[gen.dim2$row[1]:gen.dim2$row[2],gen.dim$column[1]:gen.dim$column[2]] 
  }
  
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
	gen.dim2 <- list(row=c(0,0),column=c(0,0))
	
	hp.count <- NULL; hp.mu <- NULL
	lp.count <- NULL; lp.mu <- NULL

	StopA = TRUE; StopB = TRUE
	
	ret <- list(hp.count=hp.count,lp.count=lp.count,
							hp.mu=hp.mu,lp.mu=lp.mu,
							gen.dim1=gen.dim1,gen.dim2=gen.dim2,StopA=StopA,StopB=StopB)
	ret
}

copyParentB <- function(gen.dim) {
	gen.dim1 <- list(row=c(0,0),column=c(0,0))
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
  if (part$col.part==TRUE & part$partition[2]>gen1.dim1$column[2]) STOP8
  
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
  if(part$col.part==TRUE & part$partition[2]>gen1.dim2$column[2])STOP10
  
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
												rowA1=dim1$row[1], rowA2=dim1$row[2], 
											  colA1=dim1$column[1], colA2=dim1$column[2], 
											  StopA=children$StopA, A_value=p$A_value,
											  rowB1=dim2$row[1], rowB2=dim2$row[2],
											  colB1=dim2$column[1], colB2=dim2$column[2], 
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
		row1 = curData1$rowA1[index]
		row2 = curData1$rowA2[index]
		col1 = curData1$colA1[index]
		col2 = curData1$colA2[index]
		
		A.dim <- list(row=c(row1,row2),column=c(col1,col2))

		A.counts <- counts[row1:row2, col1:col2]
		A.mu <- mu[row1:row2, col1:col2]
		
		A.Stop <- curData1$StopA[index]
		
		A_value <- curData1$A_value[index]

		row1 = curData1$rowB1[index]
		row2 = curData1$rowB2[index]
		col1 = curData1$colB1[index]
		col2 = curData1$colB2[index]
		
		B.dim <- list(row=c(row1,row2),column=c(col1,col2))
		
		B.counts <- counts[row1:row2, col1:col2]
		B.mu <- mu[row1:row2, col1:col2]
		
		B.Stop <- curData1$StopB[index]
		
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
		curData2 <- data.frame(rowA1=m$rowA1[1], A_value=formatNum(m$A_value[1]), StopA=m$StopA[1],
													 rowB1=m$rowB1[1], B_value=formatNum(m$B_value[1]), StopB=m$StopB[1])
	} else {
		index <- 0
		for (i in 1:(ite-1))
		{
			index <- index + 2^(i-1)
		}

		curData2 <- data.frame(rowA1=m$rowA1[index+1], A_value=formatNum(m$A_value[index+1]), 
													 StopA=m$StopA[index+1],
													 rowB1=m$rowB1[index+1], B_value=formatNum(m$B_value[index+1]), 
													 StopB=m$StopB[index+1])
		for (child in 2:2^(ite-1))
		{
			i <- index + child
			tempData <- data.frame(rowA1=m$rowA1[i], A_value=formatNum(m$A_value[i]), StopA=m$StopA[i],
														 rowB1=m$rowB1[i], B_value=formatNum(m$B_value[i]), StopB=m$StopB[i])
			curData2 <- rbind(curData2, tempData)
		}
	}

	return (curData2)
}

plotA <- function(x, y, p, index, colA, colP)
{
	if (p$rowA1[index] != 0) {
		if (p$StopA[index])
			text(x,y,p$A_value[index],cex=0.8,col=colP)
		else
			text(x,y,p$A_value[index],cex=0.8,col=colA)
	}
}

plotB <- function(x, y, p, index, colB, colP)
{
	if (p$rowB1[index] != 0) {
		if (p$StopB[index])
			text(x,y,p$B_value[index],cex=0.8,col=colP)
		else
			text(x,y,p$B_value[index],cex=0.8,col=colB)
	}
}

plot1Row <- function(y, p, index, colA, colB, colP)
{
	plotA(1,y,p,index,colA,colP)
	plotB(4,y,p,index,colB,colP)
	plotA(7,y,p,index+1,colA,colP)
	plotB(10,y,p,index+1,colB,colP)
}

display_values <- function(p, ite, colA, colB, colP, nrow)
{
	h <- nrow+1   # height
	
	if (ite==1) {
		y <- h*0.5
		plotA(1+1.5,y,p,1,colA,colP)
		plotB(7+1.5,y,p,1,colB,colP)
	}
	else if (ite==2) {
		y <- h*0.5; plot1Row(y,p,1,colA,colB,colP)
	}
	else if (ite==3) {
		y <- h*0.75; plot1Row(y,p,1,colA,colB,colP)
		y <- h*0.25; plot1Row(y,p,3,colA,colB,colP)
	}
	else if (ite==4) {
		y <- h*0.8; plot1Row(y,p,1,colA,colB,colP)
		y <- h*0.6; plot1Row(y,p,3,colA,colB,colP)
		y <- h*0.4; plot1Row(y,p,5,colA,colB,colP)
		y <- h*0.2; plot1Row(y,p,7,colA,colB,colP)
	}
	else if (ite==5) {
		y <- h*0.8; plot1Row(y,p,1,colA,colB,colP)
		y <- h*0.7; plot1Row(y,p,3,colA,colB,colP)
		y <- h*0.6; plot1Row(y,p,5,colA,colB,colP)
		y <- h*0.5; plot1Row(y,p,7,colA,colB,colP)

		y <- h*0.4; plot1Row(y,p,9,colA,colB,colP)
		y <- h*0.3; plot1Row(y,p,11,colA,colB,colP)
		y <- h*0.2; plot1Row(y,p,13,colA,colB,colP)
		y <- h*0.1; plot1Row(y,p,15,colA,colB,colP)
	}
	else if (ite==6) {
		y <- h*15/16; plot1Row(y,p,1,colA,colB,colP)
		y <- h*14/16; plot1Row(y,p,3,colA,colB,colP)
		y <- h*13/16; plot1Row(y,p,5,colA,colB,colP)
		y <- h*12/16; plot1Row(y,p,7,colA,colB,colP)
		
		y <- h*11/16; plot1Row(y,p,9,colA,colB,colP)
		y <- h*10/16; plot1Row(y,p,11,colA,colB,colP)
		y <- h*9/16; plot1Row(y,p,13,colA,colB,colP)
		y <- h*8/16; plot1Row(y,p,15,colA,colB,colP)

		y <- h*7/16; plot1Row(y,p,17,colA,colB,colP)
		y <- h*6/16; plot1Row(y,p,19,colA,colB,colP)
		y <- h*5/16; plot1Row(y,p,21,colA,colB,colP)
		y <- h*4/16; plot1Row(y,p,23,colA,colB,colP)
		
		y <- h*3/16; plot1Row(y,p,25,colA,colB,colP)
		y <- h*2/16; plot1Row(y,p,27,colA,colB,colP)
		y <- h*1/16; plot1Row(y,p,29,colA,colB,colP)
		y <- h*0/16; plot1Row(y,p,31,colA,colB,colP)
	}
}

plotA_withLines <- function(x, y, p, index, colA, colP, lineDown)
{
	if (p$rowA1[index] != 0) {
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

plotB_withLines <- function(x, y, p, index, colB, colP, lineDown)
{
	if (p$rowB1[index] != 0) {
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

plot2Rows <- function(y1, y2, p, index, colA, colB, colP, lineDown)
{
	plotA_withLines(2.5,y1,p,index,colA,colP,lineDown)
	plotB_withLines(8.5,y1,p,index,colB,colP,lineDown)
	plotA_withLines(2.5,y2,p,index+1,colA,colP,lineDown)
	plotB_withLines(8.5,y2,p,index+1,colB,colP,lineDown)
}

display_prev_values <- function(p, ite, colA, colB, colP, nrow)
{
	h <- nrow+1   # height
	
	if (ite==1)
	{
		y <- h*0.5+1
		plotA_withLines(2.5,y,p,1,colA,colP,0.8)
		plotB_withLines(8.5,y,p,1,colB,colP,0.8)
	} else if (ite==2)
	{
		y1 <- h*0.75+1; y2<- h*0.25+1; plot2Rows(y1,y2,p,1,colA,colB,colP,0.8)
	} else if (ite==3)
	{
		y1 <- h*0.8+1; y2 <- h*0.6+1; plot2Rows(y1,y2,p,1,colA,colB,colP,0.8)
		y1 <- h*0.4+1; y2 <- h*0.2+1; plot2Rows(y1,y2,p,3,colA,colB,colP,0.8)
	} else if (ite==4)
	{
		y1 <- h*0.8+0.6; y2 <- h*0.7+0.6; plot2Rows(y1,y2,p,1,colA,colB,colP,0.4)
		y1 <- h*0.6+0.6; y2 <- h*0.5+0.6; plot2Rows(y1,y2,p,3,colA,colB,colP,0.4)
		y1 <- h*0.4+0.6; y2 <- h*0.3+0.6; plot2Rows(y1,y2,p,5,colA,colB,colP,0.4)
		y1 <- h*0.2+0.6; y2 <- h*0.1+0.6; plot2Rows(y1,y2,p,7,colA,colB,colP,0.4)
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
	
	row1 = myData$rowA1[index]
	row2 = myData$rowA2[index]
	col1 = myData$colA1[index]
	col2 = myData$colA2[index]
	
	Stop = myData$StopA[index]
	
	rect <- list(row=c(row1,row2),column=c(col1,col2),Stop=Stop)
	
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
	
	row1 = myData$rowB1[index]
	row2 = myData$rowB2[index]
	col1 = myData$colB1[index]
	col2 = myData$colB2[index]
	
	Stop = myData$StopB[index]
	
	rect <- list(row=c(row1,row2),column=c(col1,col2),Stop=Stop)
	
	return (rect)
}

plotRectangle <- function(rect, color, nrow)
{
	m1 <- 0.04  # m1 = horizontal margin
	m2 <- 0.03  # m2 = vertical margin
	
	row1 <- rect$row[1]; row2 <- rect$row[2]
	col1 <- rect$col[1]; col2 <- rect$col[2]
	
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
