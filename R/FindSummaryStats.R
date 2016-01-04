#' @title Perform forward selection scan statistic
#' 
#' @description
#' data_xyt should contain both training data and test data.
#' Training data is mainly used for estimating means, 
#' whereas forward selection scan algorithm will be performed on the test data.
#' 
#' data_xyt should contain variables x, y and t of the same lengths.
#' (x, y) correspond to the coordinates of occurrences of incidences such as diseases. 
#' t correspond to the times of occurrences of the incidences.
#' 
#' data_poly is a data frame containing vertices of a n-sided polygon. 
#' It has n rows and 2 variables, x and y. 
#' (x, y) correspond to the coordinates of adjacent vertices of the polygon. 
#' 
#' @param data_xyt       A data frame described above.
#' @param num1Slice      Number of time units for a time slice. For instance, it is
#'                       equal to 30 for a daily time unit and a monthly time slice
#'                       (assuming 30 days per month).
#' @param num1Period     Number of time slices for a period. For instance, it is
#'                       equal to 12 for a yearly period, and a monthly time slice.
#' @param numPeriods     Number of periods required to estimate means. 
#'                       Let numTrain = num1Slice*num1Period*numPeriods. 
#'                       Then, the data with time units, t <= numTrain will be used as 
#'                       training data to estimate means.  
#'                       The data with time units, t > numTrain will be regarded as test data.
#' @param numX           Number of cells on x-axis for resulting grid.
#' @param numY           Number of cells on y-axis for resulting grid.
#' @param plotIte        Number of iterations in recursive partitioning chosen for plotting. 
#' @param cellIndices    A mx2 matrix containing (row, col) coordinates of m cells of 
#'                       resulting grid.
#'                       The time series of these m cells can be plotted later, 
#'                       using option = "specific time series" in the 
#'                       \code{\link{plotFSS}} function.
#' @param data_poly      A data frame described above.
#' @param maxSum         The algorithm would stop further partitioning 
#'                       when counts of a partition are <= maxSum.
#' @param tempSmooth     Whether to perform temporal smoothing.
#' @param spatSmooth     Whether to perform spatial smoothing. 
#' @param alpha          Parameter for temporal smoothing.
#' @param lambda         Parameter for spatial smoothing. 
#' @param num1Day        Number of time units for 1 day.
#' @param dspSlices      A vector containing time slices to be displayed later.
#' @param fillThreshold  Minimum value of z-score to highlight cells in partitions,
#'                       using option = "possible outbreaks" in the 
#'                       \code{\link{plotFSS}} function.
#' @param colPoints      Color for data points.
#' @param colA           Color for selected partitions. 
#' @param colB           Color for the other partitions.
#' @param colP           Color for pruned partitions.
#'                
#' @return Fitted model object of class "FSS". 
#'         This is assumed to be the analysis results of a function that 
#'         produces an object returned by the FSS function.
#' 
#' @export
#' 
#' @examples
#' library(RecurPart)
#' 
#' cellIndices = rbind(c(2,9),c(2,10),c(3,9),c(3,10),
#'                     c(7,5),c(7,6),c(8,6),c(8,7))
#' 
#' dspSlices <- 1:4   # a vector containing time slices to be displayed
#' 
#' # allocate the first 24 months to find monthly means of each cell, and
#' # assuming 30 days for a month
#' test <- FSS(AEGISS_ixyt, num1Slice=30, num1Period=12, numPeriods=2, numX=10, numY=9, 
#'             plotIte=8, cellIndices=cellIndices, data_poly=AEGISS_poly, lambda=0.9, 
#'             dspSlices=dspSlices, fillThreshold=0.95, colPoints="cyan")
#' 
#' plotFSS(test, option="full data", titles=c("AEGISS all data"), pause=TRUE)
#' plotFSS(test, option="overall time series", pause=TRUE)
#' plotFSS(test, option="specific time series")
#' plotFSS(test, option="possible outbreaks")
#' plotFSS(test, option="overall tree")
#' plotFSS(test, option="sub-trees")
#' plotFSS(test, option="partitions")
#' plotFSS(test, option="outbreak detection")
#' 
FSS <- function(data_xyt, num1Slice, num1Period, numPeriods, numX, numY, 
								plotIte=6, cellIndices=NULL, data_poly=NULL, 
								maxSum=3, tempSmooth=TRUE, spatSmooth=TRUE, alpha=0.7, lambda=0.8, 
								num1Day=1, dspSlices=NULL, fillThreshold=1.0, 
								colPoints="green", colA="red", colB="orange", colP="green") {
	
	numTrain <- num1Slice*num1Period*numPeriods
	
	xy <- readData(data_xyt, data_poly, numTrain)
	
	totalSlices <- ceiling(max(xy$All[,3])/num1Slice)  # total time slices for all data
	trainSlices <- num1Period*numPeriods               # num of time slices in training data
	testSlices <- totalSlices - trainSlices            # num of time slices in test data
	
	maxIte <- 6         # maxIte >= 6
	if (plotIte > maxIte) maxIte <- plotIte
	
	Obs <- convertPointsToGrid(xy, numY, numX, numTrain, num1Slice, num1Day)
	
	# checking time series of a grid cell:
	cells <- findTimeSeries(Obs, cellIndices, totalSlices, num1Period, numPeriods)
	
	outbreakData <- findOutbreakData(Obs, num1Period, numPeriods)
	
	# temporal & spatial smoothing:
	outbreakData2 <- temporalSmoothing(tempSmooth, outbreakData, alpha)
	outbreakData3 <- spatialSmoothing(spatSmooth, outbreakData2, lambda)
	
	fillColor <- "yellow"  # color used to fill rectangles in partitions, 
	                       # when z-score >= fillThreshold
	
	inputArgs <- list(numY=numY, numX=numX, 
					  				num1Slice=num1Slice, num1Period=num1Period, numPeriods=numPeriods,
						  			dspSlices=dspSlices, plotIte=plotIte,
										alpha=alpha, lambda=lambda, 
										fillThreshold=fillThreshold, fillColor=fillColor)
	
	colors <- list(colPoints=colPoints, colA=colA, colB=colB, colP=colP)
	
	partResults <- recursivePartitioning(outbreakData3, maxSum, maxIte)
	
	return (list(xy=xy, poly=data_poly, Obs=Obs, cells=cells, inputArgs=inputArgs, 
							 originalSlices=outbreakData, smoothedSlices=outbreakData3, 
							 partResults=partResults, colors=colors))
}

findOutbreakData <- function(Obs, num1Period, numPeriods) {
	nY <- dim(Obs$Counts)[1]
	nX <- dim(Obs$Counts)[2]
	nT <- dim(Obs$Counts)[3]
	
	numTrainSlices = num1Period*numPeriods         # number of slices for training data
	numTestSlices = nT - numTrainSlices            # number of slices for test data
	
	# find data and means for recursive partitioning:
	dimD = c(nY,nX,numTestSlices)
	ObsData <- array(rep(0,prod(dimD)), dimD)
	ObsMeans <- array(rep(0,prod(dimD)), dimD)
	
	means <- rep(0,num1Period)
	for (y in 1:nY)
		for (x in 1:nX) {
			ObsData[y,x,] = Obs$Counts[y,x,(numTrainSlices+1):nT]
			
			for (t in 1:num1Period) {
				sum = 0
				for (j in 1:numPeriods)
					sum = sum + Obs$Counts[y,x,t+(j-1)*num1Period]
				means[t] = sum/numPeriods
			}
			for (i in 1:numTestSlices) {
				k <- i %% num1Period
				if (k == 0) k = num1Period
				ObsMeans[y,x,i] = means[k]
			}
		}
	
	# find data and means at time 0 for temporal smoothing:
	dimD = c(nY,nX)
	ObsData0 <- array(rep(0,prod(dimD)), dimD)
	ObsMeans0 <- array(rep(0,prod(dimD)), dimD)
	
	for (y in 1:nY)
		for (x in 1:nX) {
			ObsData0[y,x] = Obs$Counts[y,x,numTrainSlices]
			ObsMeans0[y,x] = Obs$Counts[y,x,numTrainSlices]
		}
	
	return (list(data=ObsData, means=ObsMeans, data0=ObsData0, means0=ObsMeans0))
}

findTimeSeries <- function(Obs, cellIndices, totalSlices, num1Period, numPeriods) {
	numIndices <- dim(cellIndices)[1]
	
	# find sums & means:
	sums <- rep(0,numIndices)
	
	dimT = c(numIndices,totalSlices)
	means <- array(rep(0,prod(dimT)), dimT)
	timeSeries <- array(rep(0,prod(dimT)), dimT)
	
	for (i in 1:numIndices) {
		row = cellIndices[i,1]
		col = cellIndices[i,2]
		
		sum = 0
		for (t in 1:totalSlices)
			sum = sum + Obs$Counts[row,col,t]
		sums[i] = sum
		
		# find average of each month from 2 year data:
		# 		for (t in 1:12) {
		# 			means[i,t] = (Obs$Counts[row,col,t]+Obs$Counts[row,col,t+12])/2
		# 		}
		# 		for (t in 13:totalSlices) {
		# 			means[i,t] = means[i,t-12]
		# 		}
		for (t in 1:num1Period) {
			sum = 0
			for (j in 1:numPeriods)
				sum = sum + Obs$Counts[row,col,t+(j-1)*num1Period]
			means[i,t] = sum/numPeriods
		}
		for (t in (num1Period+1):totalSlices) {
			means[i,t] = means[i,t-num1Period]
		}
		
		timeSeries[i,] <- Obs$Counts[row,col,]
	}
	
	# dimensions:
	# indices: numIndices*2 (row index, col index), sums: numIndices*1, 
	# means: numIndices*totalSlices, timeSeries: numIndices*totalSlices
	return (list(indices=cellIndices, sums=sums, means=means, timeSeries=timeSeries))
}

convertPointsToGrid <- function(xy, numY, numX, numTrain, num1Slice, num1Day) {
	xmin = xy$Dim$xmin; xmax = xy$Dim$xmax; xlen = (xmax-xmin)/numX
	ymin = xy$Dim$ymin; ymax = xy$Dim$ymax; ylen = (ymax-ymin)/numY
	
	dim = c(numY,numX)   # dimension of the space (numY x numX)
	ObsAll = array(rep(0,prod(dim)), dim)     # 2D observation data
	ObsBefore = array(rep(0,prod(dim)), dim)  # 2D observation data
	ObsAfter = array(rep(0,prod(dim)), dim)   # 2D observation data
	
	npoints = dim(xy$All)[1]
	for (i in 1:npoints)
	{
		x = xy$All[i,1]
		y = xy$All[i,2]
		ix = ceiling((x-xmin)/xlen)
		iy = numY-ceiling((y-ymin)/ylen)+1
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
		iy = numY-ceiling((y-ymin)/ylen)+1
		ObsBefore[iy,ix] = ObsBefore[iy,ix]+1
	}
	
	for (i in 1:npoints)
	{
		x = xy$After[i,1]
		y = xy$After[i,2]
		if ((x <= 0) || (y <= 0)) {
			break
		}
		ix = ceiling((x-xmin)/xlen)
		iy = numY-ceiling((y-ymin)/ylen)+1
		ObsAfter[iy,ix] = ObsAfter[iy,ix]+1
	}
	
	totalSlices <- ceiling(max(xy$All[,3])/num1Slice)   # total number of time slices
	dim2 = c(numY,numX,totalSlices)                     # dimension of the space (numY x numX)
	ObsCounts = array(rep(0,prod(dim2)), dim2)          # 3D observation data
	for (i in 1:npoints)
	{
		x = xy$All[i,1]
		y = xy$All[i,2]
		t = xy$All[i,3]
		ix = ceiling((x-xmin)/xlen)
		iy = numY-ceiling((y-ymin)/ylen)+1
		it = ceiling(t/num1Slice)
		ObsCounts[iy,ix,it] = ObsCounts[iy,ix,it]+1
	}
	
	# checking counts per num1Slice time units:
	nDays = num1Day
	maxD = ceiling(max(xy$All[,3])/nDays)
	Obs1D = rep(0,maxD)
	for (i in 1:npoints)
	{
		t = xy$All[i,3]
		it = ceiling(t/nDays)
		Obs1D[it] = Obs1D[it]+1
	}
	
	nDays = num1Day*7
	maxD = ceiling(max(xy$All[,3])/nDays)
	Obs7D = rep(0,maxD)
	for (i in 1:npoints)
	{
		t = xy$All[i,3]
		it = ceiling(t/nDays)
		Obs7D[it] = Obs7D[it]+1
	}
	
	nDays = num1Day*30
	maxD = ceiling(max(xy$All[,3])/nDays)
	Obs30D = rep(0,maxD)
	for (i in 1:npoints)
	{
		t = xy$All[i,3]
		it = ceiling(t/nDays)
		Obs30D[it] = Obs30D[it]+1
	}
	
	return (list(All=ObsAll, Before=ObsBefore, After=ObsAfter, Counts=ObsCounts, 
							 Days1=Obs1D, Days7=Obs7D, Days30=Obs30D))
}

readData <- function(data_xyt, data_poly, numTrain) {
	xyAll <- cbind(data_xyt$x, data_xyt$y, data_xyt$t)
	
	xmin = min(data_xyt$x, data_poly$x); xmax = max(data_xyt$x, data_poly$x)
	ymin = min(data_xyt$y, data_poly$y); ymax = max(data_xyt$y, data_poly$y)

	m <- dim(data_xyt)[1]
	
	minT = min(data_xyt$t)
	
	vecx1 <- rep(0,m); vecy1 <- rep(0,m)    # before
	vecx2 <- rep(0,m); vecy2 <- rep(0,m)    # after
	index1 <- 1
	index2 <- 1
	for (i in 1:m)
	{
		if (data_xyt$t[i] <= (minT-1) + numTrain) {
			vecx1[index1] <- data_xyt$x[i]; vecy1[index1] <- data_xyt$y[i]
			index1 <- index1 + 1
		}
		else {
			vecx2[index2] <- data_xyt$x[i]; vecy2[index2] <- data_xyt$y[i]
			index2 <- index2 + 1
		}
	}
	xyBefore <- cbind(vecx1, vecy1)
	xyAfter <- cbind(vecx2, vecy2)
	
	# override (min, max) on x- and y-axes:
	xmin2 = findNewMin(xmin)
	xmax2 = findNewMax(xmax)
	ymin2 = findNewMin(ymin)
	ymax2 = findNewMax(ymax)
	xyDim = list(xmin=xmin2, xmax=xmax2, ymin=ymin2, ymax=ymax2)
	
	return (list(All=xyAll, Before=xyBefore, After=xyAfter, Dim=xyDim))
}

findNewMin <- function(min) {
	if (min > 0) {
		power = 10^floor(log10(min))
		newMin = floor((min/power)*10)*power/10
  }
	else {
		power = 10^floor(log10(-min))
		newMin = -floor((-min/power)*10)*power/10
	}
	return (newMin)
}

findNewMax <- function(max) {
	if (max > 0) {
		power = 10^floor(log10(max))
		newMax = ceiling((max/power)*10)*power/10
	}
	else {
		power = 10^floor(log10(-max))
		newMax = -ceiling((-max/power)*10)*power/10
	}
	return (newMax)
}
