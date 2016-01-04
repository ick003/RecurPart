### Outbreak Detection using Toy Example
#

#plot.new()   # clear plot

#############################################################################################
# defining constants for running the code:

# input data file name:
fileName = "Toy Example.csv"

mean <- 3    # expected value for each cell

# we would stop further partitioning when parent counts of a partition are <= maxSum
maxSum <- 12

# criteria = 1 (p-values), 2 (z-scores)
criteria <- 1

#lambda = 0.8   # value for spatial smoothing
lambda = 1.0   # value for spatial smoothing

# whether to plot figure:
#plotFlag <- FALSE
plotFlag <- TRUE
plotIteration <- 4       # iteration number chosen for plotting

# whether to write csv file to output p-values or z-scores:
writeCsvFlag <- TRUE

colA <- "red"     # color for selected partition
colB <- "orange"  # color for the other partition
colP <- "green"   # color for the pruned partition

#############################################################################################
# setting up data matrices:

# data matrix containing observed values
toyexample <- read.csv(fileName,sep=",",header=FALSE)
toyexample <- as.matrix(toyexample)             # the data matrix O to be partitioned

nrow = dim(toyexample)[1]
ncol = dim(toyexample)[2]

# corresponding matrix containing expected values
toyexample.mu <- matrix(rep(mean,nrow*ncol),c(nrow,ncol))

gen.dim<-list(row=c(1,nrow),column=c(1,ncol))  # dimensions of matrix to be partitioned

#############################################################################################
# spatial smoothing:

#Yt <- matrix(seq(1,2*3),c(2,3))
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
			part$col.part = FALSE; children.gen <- copyParent(retData$A.dim)
			part$values <- list(A_value=retData$A_value, B_value='')
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
			part$col.part = FALSE; children.gen <- copyParent(retData$B.dim)
			part$values <- list(A_value='', B_value=retData$B_value)
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
	par(mar=c(0,0,4,0))
	plot(c(0,11), c(0,11), type="n", xlab="", ylab="", axes=FALSE)
	
	iteStr <- ifelse(ite==1," iteration"," iterations")
	if (criteria==1)
	{
		title1 <- paste("Partitioning after ", as.character(ite), iteStr, 
										"\nwith p-values (lambda=", as.character(lambda), ")", sep="")
	} else
		title1 <- paste("Partitioning after ", as.character(ite), iteStr, 
										"\nwith z scores (lambda=", as.character(lambda), ")", sep="")
	title(main=title1, font.main=1, family="sans", cex.main=1.6)
	
	for (row in 1:10)
		for (col in 1:10)
		{
			ch <- round(toyexample[row,col])
			text(col,10-row+1,as.character(ch),cex=1)
		}
	
	# plot boundaries:
	maxColourNum = 2^ite
	colour <- rainbow(maxColourNum)
	
	for (child in 1:2^(ite-1))
	{
		rectA = getARectangle(myData, ite, child)
		# plotRectangle(rectA, colour[child])
		plotRectangle(rectA, colA)
		
		rectB = getBRectangle(myData, ite, child)
		# plotRectangle(rectB, colour[maxColourNum-child+1])
		plotRectangle(rectB, colB)
	}
	
	par(mar=c(0,0,0,0))
	plot(c(0,11), c(0,11), type="n", xlab="", ylab="", axes=FALSE)

	if (ite>1)
	{
		# display p-values of previous iteration
		values = getValues(myData, ite-1)
		display_prev_values(values, ite-1, colA, colB, colP)
	}
	
	# display p-values of current iteration
	values = getValues(myData, ite)
	display_values(values, ite, colA, colB, colP)
}







