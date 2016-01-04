
plotPossibleOutbreaks <- function(fssObj, dspSlices) {
	numY = fssObj$inputArgs$numY
	numX = fssObj$inputArgs$numX
	
	num1Slice = fssObj$inputArgs$num1Slice
	num1Period = fssObj$inputArgs$num1Period
	numPeriods = fssObj$inputArgs$numPeriods
	
	nrow = dim(fssObj$smoothedSlices$data)[1]
	ncol = dim(fssObj$smoothedSlices$data)[2]
	
	if (is.null(dspSlices)) {
		dspSlices <- fssObj$inputArgs$dspSlices
	}
	
	for (t in dspSlices)
	{
		plotOutbreakTimeSlices(fssObj$xy$All, fssObj$poly, fssObj$xy$Dim, "green", 
													 num1Slice, num1Period, numPeriods, numY, numX, 
													 fssObj$smoothedSlices$data[,,t], fssObj$smoothedSlices$means[,,t],
													 fssObj$originalSlices$data[,,t], nrow, ncol, t)

		ch = readline(prompt="Please press [Enter] to continue, or [t] to terminate plotting: ")
		if ((ch == "t") || (ch == 'T')) break
	}
}

plotOutbreakTimeSlices <- function(xyAll, poly, Dim, colP, 
																	 num1Slice, num1Period, numPeriods, numY, numX, 
																	 data, means, original, nrow, ncol, current_t) {
	myLayout <- matrix(c(1,2,3,4), nrow=2, ncol=2, byrow=TRUE)
	nf <- layout(mat=myLayout, widths=c(1,1,1,1), heights=c(1,1,1,1), respect=FALSE)
	layout.show(n = 4)
	
	xmin = Dim$xmin; xmax = Dim$xmax
	ymin = Dim$ymin;  ymax = Dim$ymax
	
	offset = num1Period*numPeriods
	
	# plot 1st graph - data:
	plotTimeSlice(xmin, xmax, ymin, ymax, xyAll, poly, num1Slice, offset, current_t, colP, 
								Dim, data, numY, numX)
	
	title <- paste("Data for time slice t=", as.character(current_t), sep="")
	title(main=title, font.main=1, family="sans", cex.main=1.6)
	
	# plot 2nd graph - means:
	plotTimeSlice(xmin, xmax, ymin, ymax, xyAll, poly, num1Slice, offset, current_t, colP, 
								Dim, means, numY, numX)
	
	title(main="Means", font.main=1, family="sans", cex.main=1.6)
	
	# plot 3rd graph - data-means>0:
	diff = data - means
	diff[diff<0] = 0
	plotTimeSlice(xmin, xmax, ymin, ymax, xyAll, poly, num1Slice, offset, current_t, colP, 
								Dim, diff, numY, numX)
	
	title(main="Data - Means > 0", font.main=1, family="sans", cex.main=1.6)
	
	# plot 4th graph - original:
	plotTimeSlice(xmin, xmax, ymin, ymax, xyAll, poly, num1Slice, offset, current_t, colP, 
								Dim, original, numY, numX)
	
	title(main="Original (before smoothing)", font.main=1, family="sans", cex.main=1.6)
}

plotTimeSlice <- function(xmin, xmax, ymin, ymax, xyAll, poly, num1Slice, offset, current_t, 
													colP, Dim, data, numY, numX) {
	par(mar=c(4,2.5,2,2))
	plot(c(xmin,xmax), c(ymin,ymax), type="n", xlab="", ylab="", axes=TRUE)
	
	plotBoundaries(poly)
	
	npoints = dim(xyAll)[1]
	for (i in 1:npoints)
	{
		t = ceiling(xyAll[i,3]/num1Slice)
		if (t == (offset+current_t)) {
			x = xyAll[i,1]
			y = xyAll[i,2]
			points(x, y, pch=16, cex=0.75, col=colP)
		}
	}
	
	plotGrid(Dim, data, numY, numX)
}

