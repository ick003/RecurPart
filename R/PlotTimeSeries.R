
# plot time series for selected cells:
plotSpecificTimeSeries <- function(fssObj, cellIndices) {
	xmin = fssObj$xy$Dim$xmin; xmax = fssObj$xy$Dim$xmax
	ymin = fssObj$xy$Dim$ymin; ymax = fssObj$xy$Dim$ymax

	myLayout <- matrix(c(1), nrow=1, ncol=1, byrow=TRUE)
	nf <- layout(mat=myLayout)

	par(mar=c(4,4,2,2))
	plot(c(xmin,xmax), c(ymin,ymax), type="n", xlab="", ylab="", axes=TRUE)
	
	num1Slice = fssObj$inputArgs$num1Slice
	num1Period = fssObj$inputArgs$num1Period
	numPeriods = fssObj$inputArgs$numPeriods
	
	if (is.null(cellIndices)) {
		cells <- fssObj$cells
	} else {
		totalSlices <- ceiling(max(fssObj$xy$All[,3])/num1Slice)
		cells <- findTimeSeries(fssObj$Obs, cellIndices, totalSlices, num1Period, numPeriods)
	}
	
	numIndices <- dim(cells$indices)[1]
	for (i in 1:numIndices) {
		cellIndex <- cells$indices[i,]
		#sum <- cells$sums[i]
		means <- cells$means[i,]
		timeSeries <- cells$timeSeries[i,]
		plotTimeSeries(cellIndex, timeSeries, means, num1Slice, num1Period, numPeriods)
		
		ch = readline(prompt="Please press [Enter] to continue, or [t] to terminate plotting: ")
		if ((ch == "t") || (ch == 'T')) break
	}
}

plotTimeSeries <- function(cellIndex, timeSeries, means, num1Slice, 
													 num1Period, numPeriods) {
	nT <- length(timeSeries)
	x = seq(1,nT)

	plot(timeSeries, pch=3, ylim=range(0,max(timeSeries)), xlab="time slice", ylab="counts")
	lines(x,timeSeries,lty=1)
	
	# plot means for each year:
	numYears = floor(nT*num1Slice/365)
	len1Year = floor(365/num1Slice)
	begin = 1; end = begin + len1Year; color = "red"
	for (i in 1:numYears) {
		points(x[begin:end],means[begin:end],pch=3,col=color)
		lines(x[begin:end],means[begin:end],lty=2,col=color)
		
		begin = end; end = begin + len1Year
		if (end > nT) end = nT
		
		# toggle color for subsequent year:
		if (color == "red") color = "blue"
		else color = "red"
	}

	# draw yearly lines:
	for (i in 1:numYears) {
		xT = ceiling(365*i/num1Slice)
		lines(c(xT,xT),c(0,max(timeSeries)),lty=2,col='green')
	}

	# plot title:
	row = cellIndex[1]
	col = cellIndex[2]
	title <- paste("Time series for (", as.character(row), ",", as.character(col), 
								 ") with sum = ", as.character(sum(timeSeries)), sep="")
	title(main=title, font.main=1, family="sans", cex.main=1.6)
	
	legend("topright", c("counts","mean"), lty=c(1,2), col=c("black","red"), 
				 text.col=c("black","red"))
	
# 	str = paste("Time series for (", as.character(row), ",", as.character(col), ")", sep="")
# 	print(str)
# 	print(timeSeries)
# 	
# 	str = paste("total sum = ",as.character(sum(timeSeries)), sep="")
# 	print(str)
	
# 	it = num1Period*numPeriods
# 	mean = 0
# 	for (t in 1:it)
# 		mean = mean + timeSeries[t]
# 	mean = mean/it
# 	str = paste("mean = ",as.character(mean), sep="")
# 	print(str)
}

plotOverallTimeSeries <- function(fssObj, titles) {
	myLayout <- matrix(c(1,2,3), nrow=3, ncol=1, byrow=TRUE)
	nf <- layout(mat=myLayout, widths=c(1,1,1), heights=c(1,1,1), respect=FALSE)
	
	plotDaySeries(fssObj$Obs$Days1, 1)
	if (is.null(titles)) title = "Time series of counts for all regions for each day"
	else title = titles[1]
	title(main=title, font.main=1, family="sans", cex.main=1.6)
	
	plotDaySeries(fssObj$Obs$Days7, 7)
	if ((is.null(titles)) || (length(titles) < 2)) title = "For each week"
	else title = titles[2]
	title(main=title, font.main=1, family="sans", cex.main=1.6)
	
	plotDaySeries(fssObj$Obs$Days30, 30)
	if ((is.null(titles)) || (length(titles) < 3)) title = "For each 30 days"
	else title = titles[3]
	title(main=title, font.main=1, family="sans", cex.main=1.6)
}

# assuming that timeSeries contains counts of nDays each
plotDaySeries <- function(timeSeries, nDays) {
	len <- length(timeSeries)
	x = seq(1,len)
	
	plot(timeSeries,pch=3,ylim=range(0,max(timeSeries)))
	lines(x,timeSeries,lty=1)
	
	numYears = ceiling(len*nDays/365)
	
	# draw yearly lines:
	for (i in 1:numYears) {
		xT = ceiling(365*i/nDays)
		lines(c(xT,xT),c(0,max(timeSeries)),lty=2,col='red')
	}
}
