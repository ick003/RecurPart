
plotOutBkDetect <- function(fssObj) {
	plotIte = fssObj$inputArgs$plotIte
	
	numY = fssObj$inputArgs$numY
	numX = fssObj$inputArgs$numX
	
	num1Slice = fssObj$inputArgs$num1Slice
	num1Period = fssObj$inputArgs$num1Period
	numPeriods = fssObj$inputArgs$numPeriods
	
	nrow = numY
	ncol = numX
	
	colA = fssObj$colors$colA
	colB = fssObj$colors$colB
	colP = fssObj$colors$colP
	
	dspSlices = fssObj$inputArgs$dspSlices
	
	alpha = fssObj$inputArgs$alpha
	lambda = fssObj$inputArgs$lambda
	
	fillThreshold = fssObj$inputArgs$fillThreshold
	fillColor = fssObj$inputArgs$fillColor
	
	for (t in dspSlices)
	{
		plotOutbreakTimeSlices2(fssObj$xy$All, fssObj$poly, fssObj$xy$Dim, 
														num1Slice, num1Period, numPeriods, numY, numX, 
														fssObj$smoothedSlices$data[,,t], fssObj$smoothedSlices$means[,,t],
														nrow, ncol, t, plotIte, fssObj$partResults, 
														colA, colB, colP, fillThreshold, fillColor)
		
		ch = readline(prompt="Please press [Enter] to continue, or [t] to terminate plotting: ")
		if ((ch == "t") || (ch == 'T')) break
	}
}

plotOutbreakTimeSlices2 <- function(xyAll, poly, Dim, 
																		num1Slice, num1Period, numPeriods, numY, numX, 
																		data, means, nrow, ncol, current_t,
																		plotIte, myData, colA, colB, colP, 
																		fillThreshold, fillColor) {
	myLayout <- matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE)
	nf <- layout(mat=myLayout, widths=c(1,1), heights=c(1,1), respect=FALSE)
	layout.show(n = 2)
	
	xmin = Dim$xmin; xmax = Dim$xmax
	ymin = Dim$ymin;  ymax = Dim$ymax
	
	offset = num1Period*numPeriods
	
	# find highlighted cells with their z scores >= fillThreshold:
	hCells <- findHighlightedCells(plotIte, nrow, ncol, myData, 
																					current_t, fillThreshold)
	
	# plot 1st graph - data-means>0:
	diff = data - means
	diff[diff<0] = 0
	plotTimeSlice2(xmin, xmax, ymin, ymax, xyAll, poly, num1Slice, offset, current_t, "green", 
								 Dim, diff, numY, numX, hCells, fillColor)
	
	title <- paste("\nData - Means > 0 \n(for time slice t=", as.character(current_t), ")", sep="")
	title(main=title, font.main=1, family="sans", cex.main=1.1)
	
	# plot 2nd graph - partition:
	dispPartitionAndValues2(plotIte, nrow, ncol, data, myData, 
												 colA, colB, colP, current_t, fillThreshold, fillColor)
}

updateHighlightedCells <- function(cells, rect, displayT, fillThreshold)
{
	row1 <- rect$nX[1]; row2 <- rect$nX[2]
	col1 <- rect$nY[1]; col2 <- rect$nY[2]
	tim1 <- rect$nT[1]; tim2 <- rect$nT[2]
	
	if ((displayT < tim1) || (tim2 < displayT))
	{
		return (cells)
	}
	
	if (row1 != 0) 
	{
		if (rect$zScore >= fillThreshold) {
			cells[row1:row2,col1:col2] = 1
		}
	}
	
	return (cells)
}

findHighlightedCells <- function(ite, nrow, ncol, myData, displayT, fillThreshold) {
	cells <- matrix(rep(0,nrow*ncol), nrow=nrow, ncol=ncol, byrow=TRUE)
	
	for (child in 1:2^(ite-1))
	{
		rectA = getARectangle(myData, ite, child)
		cells <- updateHighlightedCells(cells, rectA, displayT, fillThreshold)
		
		rectB = getBRectangle(myData, ite, child)
		cells <- updateHighlightedCells(cells, rectB, displayT, fillThreshold)
	}
	
	return (cells)
}

plotHighlightedCells <- function(Dim, numY, numX, hCells, fillColor)
{
	xmin = Dim$xmin; xmax = Dim$xmax; xlen = (xmax-xmin)/numX
	ymin = Dim$ymin; ymax = Dim$ymax; ylen = (ymax-ymin)/numY
	
	m1 <- 0.04  # m1 = horizontal margin
	m2 <- 0.03  # m2 = vertical margin
	
	# highlight cells given in hCells:
	for (row in 1:numY)
		for (col in 1:numX)
		{
			if (hCells[row, col] != 0) {
				x1 <- xmin + xlen*(col-1+m2*3); x2 <- xmin + xlen*(col-m2*3)
				y1 <- ymin + ylen*(numY-row+1-m1*1); y2 <- ymin + ylen*(numY-row+m1*1)
				
				xvec <- c(x1, x2, x2, x1)
				yvec <- c(y1, y1, y2, y2)
				polygon(xvec, yvec, col=fillColor, border = NA)
			}
		}
}

plotTimeSlice2 <- function(xmin, xmax, ymin, ymax, xyAll, poly, num1Slice, offset, current_t, 
													colP, Dim, data, numY, numX, hCells, fillColor) {
	par(mar=c(1,0,3,0))
	plot(c(xmin,xmax), c(ymin,ymax), type="n", xlab="", ylab="", axes=FALSE)
	
	plotHighlightedCells(Dim, numY, numX, hCells, fillColor)
	
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

dispPartitionAndValues2 <- function(ite, nrow, ncol, data, myData, 
																		colA, colB, colP, displayT, fillThreshold, fillColor)
{
	par(mar=c(0,0,2,0))
	plot(c(0,ncol+1), c(0,nrow+1), type="n", xlab="", ylab="", axes=FALSE)
	
	title <- paste("\n\nHighlighted cells with ",
								 "\nz scores >= ", as.character(fillThreshold), 
								 "\nafter partitioning", sep="")
	title(main=title, font.main=1, family="sans", cex.main=1.1)
	
	# plot boundaries:
	for (child in 1:2^(ite-1))
	{
		rectA = getARectangle(myData, ite, child)
		plotRectangle(rectA, colA, nrow, displayT, fillThreshold, fillColor)
		
		rectB = getBRectangle(myData, ite, child)
		plotRectangle(rectB, colB, nrow, displayT, fillThreshold, fillColor)
	}
	
	# plot numbers:
	for (row in 1:nrow)
		for (col in 1:ncol)
		{
			ch <- round(data[row,col])
			text(col,nrow-row+1,as.character(ch),cex=1)
		}
}

