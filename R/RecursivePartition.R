
recursivePartitioning <- function(outbreakData3, maxSum, maxIte) {
	toyexample <- outbreakData3$data
	toyexample.mu <- outbreakData3$means
	
	nrow = dim(toyexample)[1]
	ncol = dim(toyexample)[2]
	ntime = dim(toyexample)[3]
	
	# dimensions of array to be partitioned
	gen.dim <- list(nX=c(1,nrow), nY=c(1,ncol), nT=c(1,ntime))
	
	part <- recursive.partition3D(toyexample, toyexample.mu, gen.dim)
	children.gen <- offspring(part$col.part, part$partition, toyexample, toyexample.mu,
														gen.dim, maxSum)
	
	ite <- 1
	myData <- recordData(ite, 1, part$col.part, children.gen, part$values)
	
	curData1 <- myData
	for (ite in 2:maxIte)
	{
		curData2 <- NULL
		for (i in 1:dim(curData1)[1])
		{
			retData = getData(toyexample,toyexample.mu,curData1,i)
			
			# process A child from curData1:
			if (retData$A.Stop)
			{
				part$col.part = 1; children.gen <- copyParentA(retData$A.dim)
				#part$values <- list(A_value=retData$A_value, B_value='')
				part$values <- list(A_value=retData$A_value, B_value=as.numeric(-999))
			}
			else {
				part <- recursive.partition3D(retData$A.counts,retData$A.mu,retData$A.dim)
				children.gen <- offspring(part$col.part,part$partition,
																	toyexample,toyexample.mu,retData$A.dim,maxSum)
			}
			
			# append A partition results to curData2
			if (is.null(curData2))
				curData2 <- recordData(ite, 2*i-1, part$col.part, children.gen, part$values)
			else
			{
				tempData <- recordData(ite, 2*i-1, part$col.part, children.gen, part$values)
				curData2 <- rbind(curData2, tempData)
			}
			
			# process B child from curData1:
			if (retData$B.Stop)
			{
				part$col.part = 1; children.gen <- copyParentB(retData$B.dim)
				part$values <- list(A_value=as.numeric(-999), B_value=retData$B_value)
			}
			else {
				part <- recursive.partition3D(retData$B.counts,retData$B.mu,retData$B.dim)
				children.gen <- offspring(part$col.part,part$partition,
																	toyexample,toyexample.mu,retData$B.dim,maxSum)
			}
			
			# append B partition results to curData2
			if (is.null(curData2))
				curData2 <- recordData(ite, 2*i, part$col.part, children.gen, part$values)
			else
			{
				tempData <- recordData(ite, 2*i, part$col.part, children.gen, part$values)
				curData2 <- rbind(curData2, tempData)
			}
		}
		myData <- rbind(myData, curData2)
		curData1 <- curData2
	}

	return (myData)
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
offspring<-function(col.part, partition, scount,smu, gen.dim, maxSum) {
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
	
	return (ret)
}

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
