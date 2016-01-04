
plotOverallTree <- function(fssObj) {
	plotIte = fssObj$inputArgs$plotIte

	nrow = dim(fssObj$smoothedSlices$data)[1]
	ncol = dim(fssObj$smoothedSlices$data)[2]
	
	colA = fssObj$colors$colA
	colB = fssObj$colors$colB
	colP = fssObj$colors$colP
	
	dspSlices = fssObj$inputArgs$dspSlices
	
	if (plotIte>=6) {
		str = paste("Tree plot is not displayed for iteration > 5. ",
								"\nPlease press [Enter] to continue", sep="")
		invisible(readline(prompt=str))
		
		return (NULL)
	}
	
	if (plotIte<=0) {
		return (NULL)
	}
	
	dispWholeTree(plotIte, nrow, ncol, fssObj$partResults, colA, colB, colP, dspSlices[1], TRUE)
	invisible(readline(prompt="Please press [Enter] to continue"))
}

plotSubTrees <- function(fssObj) {
	plotIte = fssObj$inputArgs$plotIte
	
	nrow = dim(fssObj$smoothedSlices$data)[1]
	ncol = dim(fssObj$smoothedSlices$data)[2]
	
	colA = fssObj$colors$colA
	colB = fssObj$colors$colB
	colP = fssObj$colors$colP
	
	dspSlices = fssObj$inputArgs$dspSlices
	
	if (plotIte>=6) {
		str = paste("Tree plot is not displayed for iteration > 5. ",
								"\nPlease press [Enter] to continue", sep="")
		invisible(readline(prompt=str))
		
		return (NULL)
	}
	
	if (plotIte<=0) {
		return (NULL)
	}
	
	for (t in dspSlices)
	{
		dispWholeTree(plotIte, nrow, ncol, fssObj$partResults, colA, colB, colP, t, FALSE)

		ch = readline(prompt="Please press [Enter] to continue, or [t] to terminate plotting: ")
		if ((ch == "t") || (ch == 'T')) break
	}
}

dispWholeTree <- function(ite, nrow, ncol, myData, colA, colB, colP, displayT, ignoreTFlag)
{
	# configuration:
	maxIte <- ifelse(ite>5,5,ite)
	
	myLayout <- matrix(c(1), nrow=1, ncol=1, byrow=TRUE)
	nf <- layout(mat=myLayout, widths=c(1), heights=c(1), respect=FALSE)
	
	par(mar=c(0,0,1,0))
	plot(c(0,ncol+1), c(0,nrow+1), type="n", xlab="", ylab="", axes=FALSE)
	
	if (ignoreTFlag) {
		displayT <- -1
		title <- paste("z-scores of partitions for ", as.character(maxIte), 
									 " iterations for all times t", sep="")
	} else
		title <- paste("z-scores of partitions for ", as.character(maxIte), 
									 " iterations at time t=", as.character(displayT), sep="")
	title(main=title, font.main=1, family="sans", cex.main=1.1)
	
	for (i in 1:maxIte)
	{
		values = getValues(myData, i)
		displayTreeValues(values, i, colA, colB, colP, nrow, ncol, displayT, maxIte)
	}
}

getValues <- function(m, ite)
{
	if (ite==1)
	{
		curData2 <- data.frame(xA1=m$xA1[1], xA2=m$xA2[1], 
													 yA1=m$yA1[1], yA2=m$yA2[1],
													 tA1=m$tA1[1], tA2=m$tA2[1], 
													 A_value=formatNum(m$A_value[1]), StopA=m$StopA[1],
													 xB1=m$xB1[1], xB2=m$xB2[1], 
													 yB1=m$yB1[1], yB2=m$yB2[1],
													 tB1=m$tB1[1], tB2=m$tB2[1], 
													 B_value=formatNum(m$B_value[1]), StopB=m$StopB[1])
	} else {
		index <- 0
		for (i in 1:(ite-1))
		{
			index <- index + 2^(i-1)
		}
		
		curData2 <- data.frame(xA1=m$xA1[index+1], xA2=m$xA2[index+1], 
													 yA1=m$yA1[index+1], yA2=m$yA2[index+1],
													 tA1=m$tA1[index+1], tA2=m$tA2[index+1], 
													 A_value=formatNum(m$A_value[index+1]), StopA=m$StopA[index+1],
													 xB1=m$xB1[index+1], xB2=m$xB2[index+1], 
													 yB1=m$yB1[index+1], yB2=m$yB2[index+1],
													 tB1=m$tB1[index+1], tB2=m$tB2[index+1], 
													 B_value=formatNum(m$B_value[index+1]), StopB=m$StopB[index+1])
		for (child in 2:2^(ite-1))
		{
			i <- index + child
			tempData <- data.frame(xA1=m$xA1[i], xA2=m$xA2[i], 
														 yA1=m$yA1[i], yA2=m$yA2[i],
														 tA1=m$tA1[i], tA2=m$tA2[i], 
														 A_value=formatNum(m$A_value[i]), StopA=m$StopA[i],
														 xB1=m$xB1[i], xB2=m$xB2[i], 
														 yB1=m$yB1[i], yB2=m$yB2[i],
														 tB1=m$tB1[i], tB2=m$tB2[i], 
														 B_value=formatNum(m$B_value[i]), StopB=m$StopB[i])
			curData2 <- rbind(curData2, tempData)
		}
	}
	
	return (curData2)
}

formatNum <- function(n)
{
	if (n == "")
		num1 <- ""
	else
	{
		num1 <- as.numeric(n)
		num1 <- format(num1, digits=3)
	}
	return (num1)
}

displayTreeValues <- function(p, ite, colA, colB, colP, nrow, ncol, displayT, maxIte)
{
	h <- nrow+1   # height
	w <- ncol+1   # width
	
	# plot z-scores / p-values:
	x <- w*ite*0.16-w*0.18
	m <- h/2^(ite+1)
	for (i in 1:2^(ite-1))
	{
		y <- ifelse(i==1,2*m,2*m+4*m*(i-1))
		if (ite != maxIte)
			plot1HRow(x, h-y, m, p, i, colA, colB, colP, displayT)
		else
			plot1HRowWithValues(x, h-y, m, p, i, colA, colB, colP, displayT)
	}
	
	# plot joining lines for the tree:
	if (ite>1)
	{
		ox1 <- w*0.08   # offsets for line
		for (i in 1:2^(ite-1))
		{
			y <- ifelse(i==1,2*m,2*m+4*m*(i-1))
			x1 <- x-w*0.15; x2 <- x
			y2 <- h-y+m; y1 <- y2-m
			plotALine(x1+ox1,y1,x2,y2,p,i,colA,colP,displayT)
			y3 <- h-y-m
			plotBLine(x1+ox1,y1,x2,y3,p,i,colB,colP,displayT)
		}
	}
}

plot1HRow <- function(x, y, m, p, index, colA, colB, colP, displayT)
{
	plotA2(x,y+m,p,index,colA,colP,displayT)
	plotB2(x,y-m,p,index,colB,colP,displayT)
}

plotA2 <- function(x, y, p, index, colA, colP, displayT)
{
	if (displayT>0)
	{
		if ((displayT < p$tA1[index]) || (p$tA2[index] < displayT)) return (NULL)
	}
	if (p$xA1[index] != 0) {
		colorA <- ifelse(p$StopA[index],colP,colA)
		text(x,y,p$A_value[index],cex=0.8,col=colorA,pos=4)
	}
}

plotB2 <- function(x, y, p, index, colB, colP, displayT)
{
	if (displayT>0)
	{
		if ((displayT < p$tB1[index]) || (p$tB2[index] < displayT)) return (NULL)
	}
	if (p$xB1[index] != 0) {
		colorB <- ifelse(p$StopB[index],colP,colB)
		text(x,y,p$B_value[index],cex=0.8,col=colorB,pos=4)
	}
}

plot1HRowWithValues <- function(x, y, m, p, index, colA, colB, colP, displayT)
{
	plotAWithValues(x,y+m,p,index,colA,colP,displayT)
	plotBWithValues(x,y-m,p,index,colB,colP,displayT)
}

plotAWithValues <- function(x, y, p, index, colA, colP, displayT)
{
	if (y<0) return (FALSE)
	if (displayT>0)
	{
		if ((displayT < p$tA1[index]) || (p$tA2[index] < displayT)) return (FALSE)
	}
	if (p$StopA[index]) return (FALSE)
	if (p$xA1[index] != 0) {
		colorA <- ifelse(p$StopA[index],colP,colA)
		str <- paste(p$A_value[index],"  (", 
								 as.character(p$xA1[index]), ":", as.character(p$xA2[index]), ",",
								 as.character(p$yA1[index]), ":", as.character(p$yA2[index]), ")  [t=",
								 as.character(p$tA1[index]), ":", as.character(p$tA2[index]), "]", sep=" ")
		text(x,y,str,cex=0.8,col=colorA,pos=4)
		
		return (TRUE)
	} else
		return (FALSE)
}

plotBWithValues <- function(x, y, p, index, colB, colP, displayT)
{
	if (y<0) return (FALSE)
	if (displayT>0)
	{
		if ((displayT < p$tB1[index]) || (p$tB2[index] < displayT)) return (FALSE)
	}
	if (p$StopB[index]) return (FALSE)
	if (p$xB1[index] != 0) {
		colorB <- ifelse(p$StopB[index],colP,colB)
		str <- paste(p$B_value[index],"  (", 
								 as.character(p$xB1[index]), ":", as.character(p$xB2[index]), ",",
								 as.character(p$yB1[index]), ":", as.character(p$yB2[index]), ")  [t=",
								 as.character(p$tB1[index]), ":", as.character(p$tB2[index]), "]", sep=" ")
		text(x,y,str,cex=0.8,col=colorB,pos=4)
		
		return (TRUE)
	} else
		return (FALSE)
}

plotALine <- function(x1, y1, x2, y2, p, index, colA, colP, displayT)
{
	if (displayT>0)
	{
		if ((displayT < p$tA1[index]) || (p$tA2[index] < displayT)) return (NULL)
	}
	if (p$xA1[index] != 0) {
		if (p$StopA[index])
			lines(c(x1,x2),c(y1,y2),col=colP)
		else
			lines(c(x1,x2),c(y1,y2),col=colA)
	}
}

plotBLine <- function(x1, y1, x2, y2, p, index, colB, colP, displayT)
{
	if (displayT>0)
	{
		if ((displayT < p$tB1[index]) || (p$tB2[index] < displayT)) return (NULL)
	}
	if (p$xB1[index] != 0) {
		if (p$StopB[index])
			lines(c(x1,x2),c(y1,y2),col=colP)
		else
			lines(c(x1,x2),c(y1,y2),col=colB)
	}
}

