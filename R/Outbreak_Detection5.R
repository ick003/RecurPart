### Outbreak Detection: 
# 1. Generate 3D coordinates (row, col, time) directly here (without using intermediate file).
# 2. Find the 3D grid (row, col) * time.
# 3. Optionally perform temporal smoothing.
# 4. Optionally perform spatial smoothing.
# 5. Recursive (row, column and time) partitioning process on 3D grid together.
# 6. Plotting hierarchical tree of z-scores, and showing partitions for each time 1 <= t <= T.
#

#plot.new()   # clear plot

###############################################################################################
# defining constants for running the code:

X = 9
Y = 10
T = 8
displayT = 1:T      # 1 <= displayT <= endT, data at time value = displayT extracted for display
#displayT = 5:5      # 1 <= displayT <= endT, data at time value = displayT extracted for display

ignoreTFlag <- TRUE     # TRUE to ignore displayT in tree display
#ignoreTFlag <- FALSE

#plotIteration <- 9       # iteration number chosen for plotting
plotIteration <- 7       # iteration number chosen for plotting

#mean <- 3    # expected value for each cell
mean <- 3    # expected value for each cell

# define outbreak areas:
# each outbreak area (out.*) below is from rows xmin+1 to xmax, and columns from ymin+1 to ymax:
outbreak = list(nout =4, 
								out.1 = generateOutbreakList(80, 3, 7, 4, 5, 4, 5),
								out.2 = generateOutbreakList(80, 3, 7, 7, 8, 4, 5),
								out.3 = generateOutbreakList(40, 3, 4, 5, 7, 4, 5),
								out.4 = generateOutbreakList(40, 6, 7, 5, 7, 4, 5))

#############################################################################################
# the following configuration parameters seldom change:

#plotFlag <- FALSE
plotFlag <- TRUE         # whether to plot figure

temporalSmoothFlag <- TRUE
spatialSmoothFlag <- TRUE

alpha = 0.7    # value for temporal smoothing

lambda = 0.8   # value for spatial smoothing

# whether to write csv file to output p-values or z-scores:
#writeCsvFlag <- TRUE
writeCsvFlag <- FALSE

# we would stop further partitioning when parent counts of a partition are <= maxSum
maxSum <- 12*T

# criteria: only using z-scores, i.e. not using p-values
maxIte <- 6       # maxIte >= 6
if (plotIteration > maxIte) maxIte <- plotIteration

colA <- "red"     # color for selected partition
colB <- "orange"  # color for the other partition
colP <- "green"   # color for the pruned partition

#############################################################################################
# setting up 3D data array -> Obs:

# x represents rows. y represents columns:
xmin = 0; xmax = X; ymin = 0; ymax = Y; tmin = 1; tmax = T

set.seed(1)   # reinitializes the random number generator

data = generate3DPoints(mean,outbreak,xmin,xmax,ymin,ymax,tmin,tmax)

points = data$obs.all
#points = data$obs.background
#points = data$obs.outbreak

###############################################################################################
# convert 3D points to 3D grid:

Obs = convertPointsToGrid(points,xmin,xmax,ymin,ymax,tmin,tmax)

dim = c(X,Y,T)   # dimension of the space (X x Y x T)

Obs.mu <- array(rep(mean,prod(dim)), dim)

# optionally exchange column and time dimensions to test later

#############################################################################################
# temporal smoothing -> Obs2:

if (temporalSmoothFlag)
{
	Obs2 <- temporal.smooth(Obs,alpha,mean)
	
	Obs2.mu <- temporal.smooth(Obs.mu,alpha,mean)
	
	if (writeCsvFlag)
	{
		for (t in 1:T)
		{
			data = Obs2[,,t]
			filename = paste("3Ddata", as.character(t), ".csv", sep="")
			write.csv(data,file=filename,row.names=FALSE)
		}
	}
} else
{
	Obs2 <- Obs
	Obs2.mu <- Obs.mu
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
	
	if (writeCsvFlag)
	{
		for (t in 1:T)
		{
			data = Obs3[,,t]
			#data = Obs3.mu[,,t]
			filename = paste("3Ddata", as.character(t), ".csv", sep="")
			write.csv(data,file=filename,row.names=FALSE)
		}
	}
} else
{
	Obs3 <- Obs2
	Obs3.mu <- Obs2.mu
}

nrow = dim(Obs3)[1]
ncol = dim(Obs3)[2]
ntime = dim(Obs3)[3]

# dimensions of array to be partitioned
gen.dim <- list(nX=c(1,nrow), nY=c(1,ncol), nT=c(1,ntime))

#############################################################################################
# recursive partitioning process for finding outbreaks: 

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

if (writeCsvFlag)
{
	write.csv(toyexample[,,1],file="toyexample.csv",row.names=FALSE)
	write.csv(toyexample.mu[,,1],file="toyexample_mu.csv",row.names=FALSE)
	
	tempData <- myData
	cnames <- names(tempData)
	cnames[11] <- 'A z-score'; cnames[19] <- 'B z-score'
	names(tempData) <- cnames
	write.csv(tempData,file="z scores.csv",row.names=FALSE)
}

if (plotFlag)
{
	# plot points and counts:
	colBoundary <- "orange"; 	colText <- "red"
	plot3DPoints(points,Obs,xmin,xmax,ymin,ymax,tmin,tmax,colBoundary,colText)
	
	# plot hierarchical trees:
	plot.new()   # clear plot
	
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







