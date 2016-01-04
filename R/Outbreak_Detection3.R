### Outbreak Detection using 3D data
# (row and column partitioning on 2D only, i.e. partitioning independent of other times)
#

#plot.new()   # clear plot

#############################################################################################
# defining constants for running the code:

# input data file name:
fileName = "3DdataTmp.csv"

X = 9
Y = 10
T = 8

mean <- 3    # expected value for each cell

# we would stop further partitioning when parent counts of a partition are <= maxSum
maxSum <- 12

# criteria = 1 (p-values), 2 (z-scores)
criteria <- 2

alpha = 0.7    # value for temporal smoothing

lambda = 0.8   # value for spatial smoothing

partT = 1      # 1 <= partT <= T, data at time value = partT extracted for partitioning

# whether to plot figure:
#plotFlag <- FALSE
plotFlag <- TRUE
plotIteration <- 5       # iteration number chosen for plotting

# whether to write csv file to output p-values or z-scores:
writeCsvFlag <- TRUE

colA <- "red"     # color for selected partition
colB <- "orange"  # color for the other partition
colP <- "green"   # color for the pruned partition

#############################################################################################
# setting up 3D data array:

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

#############################################################################################
# temporal smoothing:

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

# data matrix containing observed values
toyexample <- Obs2[,,partT]
toyexample <- as.matrix(toyexample)             # the data matrix to be partitioned

nrow = dim(toyexample)[1]
ncol = dim(toyexample)[2]

# corresponding matrix containing expected values
toyexample.mu <- Obs2.mu[,,partT]
toyexample.mu <- as.matrix(toyexample.mu)

gen.dim<-list(row=c(1,nrow),column=c(1,ncol))  # dimensions of matrix to be partitioned

#############################################################################################
# spatial smoothing:

toyexample <- spatial.smooth(toyexample,lambda)

toyexample.mu <- spatial.smooth(toyexample.mu,lambda)

#############################################################################################
# recursive partitioning process for finding outbreaks: 

part <- recursive.partition(criteria,toyexample,toyexample.mu,gen.dim)
children.gen <- offspring(part$col.part,part$partition,toyexample,toyexample.mu,gen.dim,maxSum)

ite <- 1
myData <- recordData(ite, 1, part$col.part, children.gen, part$values)

curData1 <- myData
for (ite in 2:6)
{
	curData2 <- NULL
	for (i in 1:dim(curData1)[1])
	{
		retData = getData(toyexample,toyexample.mu,curData1,i)

		# process A child from curData1:
		if (retData$A.Stop)
		{
			part$col.part = FALSE; children.gen <- copyParentA(retData$A.dim)
			#part$values <- list(A_value=retData$A_value, B_value='')
			part$values <- list(A_value=retData$A_value, B_value=as.numeric(-999))
		}
		else {
			part <- recursive.partition(criteria,retData$A.counts,retData$A.mu,retData$A.dim)
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
			part$col.part = FALSE; children.gen <- copyParentB(retData$B.dim)
			part$values <- list(A_value=as.numeric(-999), B_value=retData$B_value)
		}
		else {
			part <- recursive.partition(criteria,retData$B.counts,retData$B.mu,retData$B.dim)
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
	write.csv(toyexample,file="toyexample.csv",row.names=FALSE)
	write.csv(toyexample.mu,file="toyexample_mu.csv",row.names=FALSE)
	
	tempData <- myData
	cnames <- names(tempData)
	if (criteria==1) {
		cnames[9] <- 'A p-value'; cnames[15] <- 'B p-value'
		names(tempData) <- cnames
		write.csv(tempData,file="p-values.csv",row.names=FALSE)
	} else {
		cnames[9] <- 'A z-score'; cnames[15] <- 'B z-score'
		names(tempData) <- cnames
		write.csv(tempData,file="z scores.csv",row.names=FALSE)
	}
}

if (plotFlag)
{
	ite <- plotIteration
	
	myLayout <- matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE)
	nf <- layout(mat=myLayout, widths=c(4,3), heights=c(1,1), respect=FALSE)
	
	# plot numbers:
	par(mar=c(0,0,5,0))
	plot(c(0,ncol+1), c(0,nrow+1), type="n", xlab="", ylab="", axes=FALSE)
	
	iteStr <- ifelse(ite==1," iteration"," iterations")
	criStr <- ifelse(criteria==1, "p-values", "z scores")
	title1 <- paste("Partitioning after ", as.character(ite), iteStr, 
									"\nwith ", criStr, " (t=", as.character(partT), ",", 
									"\nalpha=", as.character(alpha), ", ", 
									"lambda=", as.character(lambda), ")", sep="")
	title(main=title1, font.main=1, family="sans", cex.main=1.6)
	
	for (row in 1:nrow)
		for (col in 1:ncol)
		{
			ch <- round(toyexample[row,col])
			text(col,nrow-row+1,as.character(ch),cex=1)
		}
	
	# plot boundaries:
	maxColourNum = 2^ite
	colour <- rainbow(maxColourNum)
	
	for (child in 1:2^(ite-1))
	{
		rectA = getARectangle(myData, ite, child)
		# plotRectangle(rectA, colour[child])
		plotRectangle(rectA, colA, nrow)
		
		rectB = getBRectangle(myData, ite, child)
		# plotRectangle(rectB, colour[maxColourNum-child+1])
		plotRectangle(rectB, colB, nrow)
	}
	
	par(mar=c(0,0,0,0))
	plot(c(0,ncol+1), c(0,nrow+1), type="n", xlab="", ylab="", axes=FALSE)

	if (ite>1)
	{
		# display p-values of previous iteration
		values = getValues(myData, ite-1)
		display_prev_values(values, ite-1, colA, colB, colP, nrow)
	}
	
	# display p-values of current iteration
	values = getValues(myData, ite)
	display_values(values, ite, colA, colB, colP, nrow)
}







