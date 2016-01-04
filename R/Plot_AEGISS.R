# plot AEGISS dot pattern for time: (T1 <= t <= T2) and (t > T2) respectively
#
# AEGISS_ixyt.txt = patient id, x coordinate, y coordinate, time t
#   10572 cases of non-specific gastrointestinal disease in the country of Hampshire, UK,
#     as reported to NHS Direct.
#   Each (x,y)-location corresponds to the centroid of the unit post-code of the residental address
#     of the person making the call to NHS Direct. Unit of distance = 1 km.
#   Unit of time is 1 day, with day 1 corresponding to 1 Jan, 2001.
#
# AEGISS_poly.txt = xy coordinates of 120-sided polygon representing the boundary of the study-region
#   (Hampshire, UK). Unit of distance = 1 km.
#
#

# specifying dividing time to plot (T1 <= t <= T2) and (t > T2) respectively
T1 = 1
T2 = 30*24     # allocate the first 24 months to find monthly means of each cell

xy <- readAegissData("data\\AEGISS_ixyt.txt", T1, T2)

nD = 30        # number of days for each observation, i.e. for a time slice
#nD = 14        # number of days for each observation, i.e. for a time slice
#nD = 7        # number of days for each observation, i.e. for a time slice

nT <- ceiling(max(xy$All[,3])/nD)          # number of time slices for all 3 years = 37
T <- nT - ceiling(T2/nD)                   # number of time slices for data to detect outbreaks = 13

#fine = TRUE
fine = FALSE

if (fine)
{
	nY = 20      # number of divisions on y-axis
	nX = 20      # number of divisions on x-axis
} else {
	nY = 9       # number of divisions on y-axis
	nX = 10      # number of divisions on x-axis
}

# cell indices (row, col) used to find and plot time series:
# cellIndices = rbind(c(6,5),c(6,6),c(6,7),c(6,9),
# 										c(7,4),c(7,5),c(7,6),c(7,7),c(7,8),
# 										c(1,7),
# 										c(2,7),c(2,8),c(2,9),c(2,10),
# 										c(3,4),c(3,5),c(3,7),c(3,10),
# 										c(4,9),c(4,10),
# 										c(5,6),c(5,9),
# 										c(8,5),c(8,6),c(8,7),c(8,8),c(9,3),c(9,4))
cellIndices = rbind(c(2,9),c(2,10),c(3,9),c(3,10),
										c(7,5),c(7,6),c(8,6),c(8,7))
#cellIndices = rbind(c(6,14))

colPoints <- "green"   # color for the data points

###################################################################################################

#displayT = 1:T      # 1 <= displayT <= endT, data at time value = displayT extracted for display
#displayT = 5:5      # 1 <= displayT <= endT, data at time value = displayT extracted for display
#displayT = 5:6
displayT = 1:4

ignoreTFlag <- TRUE     # TRUE to ignore displayT in tree display
#ignoreTFlag <- FALSE

#plotTimeSeriesFlag <- FALSE
plotTimeSeriesFlag <- TRUE

#analysisFlag <- FALSE
analysisFlag <- TRUE

#plotFlag <- FALSE
plotFlag <- TRUE         # whether to plot figure

plotIteration <- 8       # iteration number chosen for plotting

# we would stop further partitioning when parent counts of a partition are <= maxSum
#maxSum <- 12*T
maxSum <- 3

# criteria: only using z-scores, i.e. not using p-values
maxIte <- 6       # maxIte >= 6
if (plotIteration > maxIte) maxIte <- plotIteration

temporalSmoothFlag <- TRUE
spatialSmoothFlag <- TRUE

alpha = 0.7    # value for temporal smoothing

#lambda = 0.8   # value for spatial smoothing
lambda = 0.9   # value for spatial smoothing

colA <- "red"     # color for selected partition
colB <- "orange"  # color for the other patrtition
colP <- "green"   # color for the pruned partition

###################################################################################################

#library(splancs)

Obs <- convertAegissPointsToGrid(xy, nY, nX, T2, nD)

# checking time series of a grid cell:
all <- Obs$All         # each cell in (nY*nX) grid contains total counts for all 3 years
cells <- findTimeSeries(Obs,cellIndices,nT,nD)

outbreakData <- findOutbreakData(Obs, T2, nD)
means2D <- outbreakData$means[,,1]

poly <- readAegissPoly("data\\AEGISS_poly.txt")

plot.new()

# plotting time series of grid cells:
if (plotTimeSeriesFlag) {
	# plotting poly:
	title = "AEGISS all data"
	plotAegiss2(xy$All, poly, xy$Dim, title, Obs, nY, nX, colPoints, TRUE)
	#plotAegiss2(xy$All, poly, xy$Dim, title, Obs, nY, nX, colPoints, FALSE)

	# plot time series of all regions for each day:
	plotDaySeries(Obs)
	
	ch = readline(prompt="Press [enter] to continue, or [t] to terminate: ")
	myLayout <- matrix(c(1,1), nrow=1, ncol=1, byrow=TRUE)
	nf <- layout(mat=myLayout)

	if ((ch != "t") && (ch != 'T')) {
		# plot time series for selected cells:
		numIndices <- dim(cellIndices)[1]
		for (i in 1:numIndices) {
			cellIndex <- cells$indices[i,]
			sum <- cells$sums[i]
			means <- cells$means[i,]
			timeSeries <- cells$timeSeries[i,]
			plotTimeSeries(timeSeries,T2,nD,cellIndex,means)
			ch = readline(prompt="Press [enter] to continue, or [t] to terminate: ")
			if ((ch == "t") || (ch == 'T')) break
		}
	}
	else
		analysisFlag = FALSE
}

# # plotting all points:
# title = "AEGISS all data"
# plotAegiss(xy$All, poly, title, TRUE)
# 
# # plotting points <= time T:
# title = paste("AEGISS data for ", as.character(T1), "<= t <= ", as.character(T2), sep="")
# plotAegiss(xy$Before, poly, title, TRUE)
# 
# # plotting points > time T:
# title = paste("AEGISS data for t > ", as.character(T2), sep="")
# plotAegiss(xy$After, poly, title, FALSE)

#############################################################################################

dim = dim(outbreakData$data)   # dimension of the space (Y x X x T)
T = dim[3]

# temporal smoothing -> Obs2:

if (temporalSmoothFlag)
{
	Obs2 <- temporal.smooth2(outbreakData$data,alpha,outbreakData$data0)
	
	Obs2.mu <- temporal.smooth2(outbreakData$means,alpha,outbreakData$means0)
} else {
	Obs2 <- outbreakData$data
	Obs2.mu <- outbreakData$means
}

#############################################################################################
# spatial smoothing -> Obs3:

if (spatialSmoothFlag)
{
	Obs3 <- array(rep(0,prod(dim)), dim)
	Obs3.mu <- array(rep(0,prod(dim)), dim)
	
	for (t in 1:T)
	{
		# data matrix containing observed values
		data <- Obs2[,,t]
		data <- as.matrix(data)             # the data matrix to be partitioned
		Obs3[,,t] <- spatial.smooth(data,lambda)
		
		# corresponding matrix containing expected values
		data.mu <- Obs2.mu[,,t]
		data.mu <- as.matrix(data.mu)
		Obs3.mu[,,t] <- spatial.smooth(data.mu,lambda)
	}
} else {
	Obs3 <- Obs2
	Obs3.mu <- Obs2.mu
}

nrow = dim(Obs3)[1]
ncol = dim(Obs3)[2]
ntime = dim(Obs3)[3]

# dimensions of array to be partitioned
gen.dim <- list(nX=c(1,nrow), nY=c(1,ncol), nT=c(1,ntime))

#############################################################################################

# plot data + (data-means) in grid for each time slice:
# Obs$Counts[,,24+t] = outbreakData$data[,,t] vs outbreakData$means[,,t]
# Obs2[,,t] vs Obs2.mu[,,t]
# Obs3[,,t] vs Obs3.mu[,,t]
for (t in displayT)
{
	plotOutbreakTimeSlices(xy$All, poly, xy$Dim, "green", T2, nD, nY, nX, 
												 Obs3[,,t], Obs3.mu[,,t], outbreakData$data[,,t], nrow, ncol, t)
	ch = readline(prompt="Press [enter] to continue, or [t] to terminate: ")
	if ((ch == "t") || (ch == 'T')) {
		analysisFlag <- FALSE
		break
	}
}

#############################################################################################
# recursive partitioning process for finding outbreaks: 

if (analysisFlag) {
	toyexample <- Obs3
	toyexample.mu <- Obs3.mu
	
	part <- recursive.partition3D(toyexample,toyexample.mu,gen.dim)
	children.gen <- offspring(part$col.part,part$partition,toyexample,toyexample.mu,gen.dim,maxSum)
	
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
}

if (analysisFlag && plotFlag)
{
	# plot hierarchical trees:
	#plot.new()   # clear plot
	
	if (plotIteration<6)
	{
		dispWholeTree(plotIteration, nrow, ncol, myData, colA, colB, colP, displayT[1], TRUE)
		invisible(readline(prompt="Press [enter] to continue"))
		if (!ignoreTFlag) {
			for (t in displayT)
			{
				dispWholeTree(plotIteration, nrow, ncol, myData, colA, colB, colP, t, ignoreTFlag)
				invisible(readline(prompt="Press [enter] to continue"))
			}
		}
	} else {
		str = paste("Tree plot is not displayed for iteration > 5. ",
								"\nPlease press [enter] to continue", sep="")
		invisible(readline(prompt=str))
	}
	
	for (t in displayT)
	{
		dispPartitionAndValues(plotIteration, nrow, ncol, alpha, lambda, toyexample, myData, 
													 colA, colB, colP, t)
		if (t<T) {
			ch = readline(prompt="Press [enter] to continue, or [t] to terminate: ")
			if ((ch == "t") || (ch == 'T')) break
		}
	}
	
	#dispPartitionAndTree(plotIteration, nrow, ncol, alpha, lambda, toyexample, myData, 
	#										 colA, colB, colP, displayT)
}
