
plotPartitions <- function(fssObj) {
	plotIte = fssObj$inputArgs$plotIte
	
	nrow = dim(fssObj$smoothedSlices$data)[1]
	ncol = dim(fssObj$smoothedSlices$data)[2]
	
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
		dispPartitionAndValues(plotIte, nrow, ncol, alpha, lambda, 
													 fssObj$smoothedSlices$data, fssObj$partResults, 
													 colA, colB, colP, t, fillThreshold, fillColor)

		ch = readline(prompt="Please press [Enter] to continue, or [t] to terminate plotting: ")
		if ((ch == "t") || (ch == 'T')) break
	}
}

dispPartitionAndValues <- function(ite, nrow, ncol, alpha, lambda, toyexample, myData, 
																	 colA, colB, colP, displayT, fillThreshold, fillColor)
{
	myLayout <- matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE)
	nf <- layout(mat=myLayout, widths=c(4,3), heights=c(1,1), respect=FALSE)
	
	par(mar=c(0,0,5,0))
	plot(c(0,ncol+1), c(0,nrow+1), type="n", xlab="", ylab="", axes=FALSE)
	
	iteStr <- ifelse(ite==1," iteration"," iterations")
	title <- paste("Partitioning after ", as.character(ite), iteStr, 
								 "\nwith z scores (t=", as.character(displayT), ",", 
								 "\nalpha=", as.character(alpha), ", ", 
								 "lambda=", as.character(lambda), ")", sep="")
	title(main=title, font.main=1, family="sans", cex.main=1.6)
	
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
			ch <- round(toyexample[row,col,displayT])
			text(col,nrow-row+1,as.character(ch),cex=1)
		}
	
	par(mar=c(0,0,5,0))
	plot(c(0,ncol+1), c(0,nrow+1), type="n", xlab="", ylab="", axes=FALSE)
	
	title <- paste("z-score table",
								 "\n(row1:row2,col1:col2) ", 
								 "\n[t=tim1:tim2]", sep="")
	title(main=title, font.main=1, family="sans", cex.main=1.6)
	
	# display z-scores of current iteration
	values = getValues2(myData, ite)
	sortedValues = sortZScores(values, ite, colA, colB, colP, displayT)
	displayValuesInRows(sortedValues, ite, colA, colB, colP, nrow, ncol, displayT)
}

getValues2 <- function(m, ite)
{
	if (ite==1)
	{
		curData2 <- data.frame(xA1=m$xA1[1], xA2=m$xA2[1], 
													 yA1=m$yA1[1], yA2=m$yA2[1],
													 tA1=m$tA1[1], tA2=m$tA2[1], 
													 A_value=m$A_value[1], StopA=m$StopA[1],
													 xB1=m$xB1[1], xB2=m$xB2[1], 
													 yB1=m$yB1[1], yB2=m$yB2[1],
													 tB1=m$tB1[1], tB2=m$tB2[1], 
													 B_value=m$B_value[1], StopB=m$StopB[1])
	} else {
		index <- 0
		for (i in 1:(ite-1))
		{
			index <- index + 2^(i-1)
		}
		
		curData2 <- data.frame(xA1=m$xA1[index+1], xA2=m$xA2[index+1], 
													 yA1=m$yA1[index+1], yA2=m$yA2[index+1],
													 tA1=m$tA1[index+1], tA2=m$tA2[index+1], 
													 A_value=m$A_value[index+1], StopA=m$StopA[index+1],
													 xB1=m$xB1[index+1], xB2=m$xB2[index+1], 
													 yB1=m$yB1[index+1], yB2=m$yB2[index+1],
													 tB1=m$tB1[index+1], tB2=m$tB2[index+1], 
													 B_value=m$B_value[index+1], StopB=m$StopB[index+1])
		for (child in 2:2^(ite-1))
		{
			i <- index + child
			tempData <- data.frame(xA1=m$xA1[i], xA2=m$xA2[i], 
														 yA1=m$yA1[i], yA2=m$yA2[i],
														 tA1=m$tA1[i], tA2=m$tA2[i], 
														 A_value=m$A_value[i], StopA=m$StopA[i],
														 xB1=m$xB1[i], xB2=m$xB2[i], 
														 yB1=m$yB1[i], yB2=m$yB2[i],
														 tB1=m$tB1[i], tB2=m$tB2[i], 
														 B_value=m$B_value[i], StopB=m$StopB[i])
			curData2 <- rbind(curData2, tempData)
		}
	}
	
	return (curData2)
}

getARectangle <- function(myData, ite, child)
{
	index <- 0
	if (ite>1)
	{
		for (i in 1:(ite-1))
		{
			index <- index + 2^(i-1)
		}
	}
	index <- index + child
	
	row1 = myData$xA1[index]
	row2 = myData$xA2[index]
	col1 = myData$yA1[index]
	col2 = myData$yA2[index]
	tim1 = myData$tA1[index]
	tim2 = myData$tA2[index]
	
	Stop = myData$StopA[index]
	
	zScore = myData$A_value[index]
	
	rect <- list(nX=c(row1,row2),nY=c(col1,col2),nT=c(tim1,tim2),Stop=Stop,zScore=zScore)
	
	return (rect)
}

getBRectangle <- function(myData, ite, child)
{
	index = 0
	if (ite>1)
	{
		for (i in 1:(ite-1))
		{
			index <- index + 2^(i-1)
		}
	}
	index <- index + child
	
	row1 = myData$xB1[index]
	row2 = myData$xB2[index]
	col1 = myData$yB1[index]
	col2 = myData$yB2[index]
	tim1 = myData$tB1[index]
	tim2 = myData$tB2[index]
	
	Stop = myData$StopB[index]
	
	zScore = myData$B_value[index]
	
	rect <- list(nX=c(row1,row2),nY=c(col1,col2),nT=c(tim1,tim2),Stop=Stop,zScore=zScore)
	
	return (rect)
}

plotRectangle <- function(rect, color, nrow, displayT, fillThreshold, fillColor)
{
	m1 <- 0.04  # m1 = horizontal margin
	m2 <- 0.03  # m2 = vertical margin
	
	row1 <- rect$nX[1]; row2 <- rect$nX[2]
	col1 <- rect$nY[1]; col2 <- rect$nY[2]
	tim1 <- rect$nT[1]; tim2 <- rect$nT[2]
	
	if ((displayT < tim1) || (tim2 < displayT))
	{
		return ()
	}
	
	if (rect$Stop) color <- "green"
	
	if (row1 != 0) 
	{
		x1 <- col1-0.5+m2; x2 <- col2+0.5-m2
		y1 <- nrow-row1+1+0.5-m1; y2 <- nrow-row2+1-0.5+m1
		
		# top line:
		x <- c(x1, x2); y <- c(y1,y1)
		lines(x,y,col=color)
		
		# bottom line:
		x <- c(x1, x2); y <- c(y2,y2)
		lines(x,y,col=color)
		
		# left line:
		x <- c(x1, x1); y <- c(y1,y2)
		lines(x,y,col=color)
		
		# right line:
		x <- c(x2, x2); y <- c(y1,y2)
		lines(x,y,col=color)
		
		if (rect$zScore >= fillThreshold) {
			x1 <- x1+m2*3; x2 <- x2-m2*3
			y1 <- y1-m1*2; y2 <- y2+m1*2
			
			xvec <- c(x1, x2, x2, x1)
			yvec <- c(y1, y1, y2, y2)
			polygon(xvec, yvec, col=fillColor, border = NA)
		}
	}
}

sortZScores <- function(p, ite, colA, colB, colP, displayT)
{
	# sort z-scores (rows:cols) [t=times]:
	zScoreData <- NULL
	for (i in 1:2^(ite-1))
	{
		curZScore <- getA_ZScore(p,i,colA,colP,displayT)
		zScoreData <- rbind(zScoreData, curZScore)
		curZScore <- getB_ZScore(p,i,colB,colP,displayT)
		zScoreData <- rbind(zScoreData, curZScore)
	}
	
	# sort by 1st column in decreasing order:
	sorted <- zScoreData[ order(-zScoreData[,1]), ]
	
	return (sorted)
}

getA_ZScore <- function(p, index, colA, colP, displayT)
{
	if (displayT>0)
	{
		if ((displayT < p$tA1[index]) || (p$tA2[index] < displayT)) return (NULL)
	}
	if (p$StopA[index]) return (NULL)
	if (p$xA1[index] != 0) {
		colorA <- ifelse(p$StopA[index],colP,colA)

		return (data.frame(zscore=p$A_value[index], 
											 x1=p$xA1[index], x2=p$xA2[index], 
											 y1=p$yA1[index], y2=p$yA2[index], 
											 t1=p$tA1[index], t2=p$tA2[index], color=colorA))
	} else
		return (NULL)
}

getB_ZScore <- function(p, index, colB, colP, displayT)
{
	if (displayT>0)
	{
		if ((displayT < p$tB1[index]) || (p$tB2[index] < displayT)) return (NULL)
	}
	if (p$StopB[index]) return (NULL)
	if (p$xB1[index] != 0) {
		colorB <- ifelse(p$StopB[index],colP,colB)
		
		return (data.frame(zscore=p$B_value[index], 
											 x1=p$xB1[index], x2=p$xB2[index], 
											 y1=p$yB1[index], y2=p$yB2[index], 
											 t1=p$tB1[index], t2=p$tB2[index], color=colorB))
	} else
		return (NULL)
}

displayValuesInRows <- function(p, ite, colA, colB, colP, nrow, ncol, displayT)
{
	h <- nrow+1   # height
	w <- ncol+1   # width
	
	# plot z-scores (rows:cols) [t=times]:
	x <- w*0.01
	m <- h/35
	y <- 2*m
	for (i in 1:dim(p)[1])
	{
		flag <- plotZScores(x,h-y,p,i)
		if (flag) y <- y+2*m
		else break
	}
}

plotZScores <- function(x, y, p, index)
{
	if (y<0) return (FALSE)
	
	z <- format(p$zscore[index], digits=3)
	
	str <- paste(z,"  (", 
							 as.character(p$x1[index]), ":", as.character(p$x2[index]), ",",
							 as.character(p$y1[index]), ":", as.character(p$y2[index]), ")  [t=",
							 as.character(p$t1[index]), ":", as.character(p$t2[index]), "]", sep=" ")
	text(x,y,str,cex=0.8,col=as.character(p$color[index]),pos=4)
	
	return (TRUE)
}

