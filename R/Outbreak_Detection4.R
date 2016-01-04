### Outbreak Detection using 3D data
# (row, column and time partitioning on 3D together)
#

#plot.new()   # clear plot

#############################################################################################
# defining constants for running the code:

# input data file name:

fileName = "3DdataTmp.csv"
X = 9
Y = 10
T = 8
displayT = 1:T      # 1 <= displayT <= endT, data at time value = displayT extracted for display
#displayT = 5:5      # 1 <= displayT <= endT, data at time value = displayT extracted for display
temporalSmoothFlag <- TRUE
spatialSmoothFlag <- TRUE

ignoreTFlag <- TRUE     # TRUE to ignore displayT in tree display
#ignoreTFlag <- FALSE

# whether to plot figure:
#plotFlag <- FALSE
plotFlag <- TRUE
plotIteration <- 14       # iteration number chosen for plotting

# fileName = "3DdataTmp2.csv"
# X = 3
# Y = 3
# T = 3
# displayT = 3      # 1 <= displayT <= T, data at time value = displayT extracted for display
# temporalSmoothFlag <- FALSE
# spatialSmoothFlag <- FALSE

# fileName = "Toy Example 3D.csv"
# X = 10
# Y = 10
# T = 6
# displayT = 6      # 1 <= displayT <= T, data at time value = displayT extracted for display
# temporalSmoothFlag <- FALSE
# spatialSmoothFlag <- FALSE

mean <- 3    # expected value for each cell

# we would stop further partitioning when parent counts of a partition are <= maxSum
maxSum <- 12*T

# criteria: only using z-scores, i.e. not using p-values

alpha = 0.7    # value for temporal smoothing

lambda = 0.8   # value for spatial smoothing

# whether to write csv file to output p-values or z-scores:
#writeCsvFlag <- TRUE
writeCsvFlag <- FALSE

colA <- "red"     # color for selected partition
colB <- "orange"  # color for the other partition
colP <- "green"   # color for the pruned partition

maxIte <- 14       # maxIte >= 6

#############################################################################################
# setting up 3D data array -> Obs:

dim = c(X,Y,T)   # dimension of the space (X x Y x T)

ObsTmp <- read.csv(fileName,sep=",",header=TRUE)
ObsTmp <- as.matrix(ObsTmp)

Obs <- array(rep(0,prod(dim)), dim)
for (i in 1:T)
{
	col1 = (i-1)*Y+1
	col2 = (i-1)*Y+Y
	Obs[1:X, 1:Y, i] = ObsTmp[1:X, col1:col2]
}

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
		invisible(readline(prompt="Press [enter] to continue"))
	}
											 
 	#dispPartitionAndTree(plotIteration, nrow, ncol, alpha, lambda, toyexample, myData, 
 	#										 colA, colB, colP, displayT)
}







