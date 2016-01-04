#' @title Plot FSS analysis results
#' 
#' @description
#' Plot analysis results after calling the \code{\link{FSS}} function. 
#' 
#' The "option" parameter must be one of the following values:
#' 
#' "full data" - plotting data points and the n-sided polygon (if available), 
#' and the corresponding cell counts in grid.
#' 
#' "overall time series" - plotting daily, weekly and monthly 
#' (assuming 30 days per month) time series of all the counts 
#' in all the regions so that it may be easier to observe 
#' any periodic or seasonal trend in the data.
#' 
#' "specific time series" - plotting time series of counts for selected cells, 
#' given by cellIndices. 
#' Besides, the corresponding expected means will also be displayed.
#' The "titles" and pause" parameters are ignored when this option is selected.
#' 
#' "possible outbreaks" - plotting the following 4 figures per time slice:
#' Fig 1 (top left)     - data after temporal and spatial smoothing.
#' Fig 2 (top right)    - expected counts.
#' Fig 3 (bottom left)  - the positive difference between counts in Fig 1 and Fig 2.
#' Fig 4 (bottom right) - original data before temporal and spatial smoothing. 
#' The time slices are given by "dspSlices" in the \code{\link{FSS}} function. 
#' The "titles" and pause" parameters are ignored when this option is selected.
#' 
#' "overall trees" - plotting a hierarchical tree of z-scores for all time slices.
#' The number of iterations is given by "plotIte" in the \code{\link{FSS}} function. 
#' "plotIte" must be less than 6, in order to plot this graph. 
#' The "titles" and pause" parameters are ignored when this option is selected.
#' 
#' "sub-trees" - plotting hierarchical trees of z-scores for individual time slices,
#' given by "dspSlices". The number of iterations is given by "plotIte". 
#' Both "dspSlices" and "plotIte" can be specified in the \code{\link{FSS}} function. 
#' "plotIte" must be less than 6, in order to plot this graph. 
#' The "titles" and pause" parameters are ignored when this option is selected.
#' 
#' "partitions" - plotting the results of partitions and their z-scores for 
#' individual time slices. The highlighted yellow partitions represent possible
#' outbreaks, with z-scores >= "fillThreshold". 
#' The time slice numbers are given by "dspSlices". 
#' The number of iterations to perform partitioning is given by "plotIte". 
#' "fillThreshold", "dspSlices" and "plotIte" can be specified in the 
#' \code{\link{FSS}} function. 
#' The "titles" and pause" parameters are ignored when this option is selected.
#' 
#' "outbreak detection" - plotting to highlight areas of possible outbreaks. 
#' It will plot the following 2 figures per time slice:
#' Fig 1 (left)  - Data - Mean > 0. It is similar to Fig.3 for 
#' "option" = "possible outbreaks". 
#' Fig 2 (right) - It is similar to the figure for "option" = "partitions". 
#' The cells with z-score greater than "fillThreshold" are highlighted both in 
#' Fig 1 and Fig 2.
#' "fillThreshold" and "dspSlices" can be specified in the \code{\link{FSS}} function. 
#' The "titles" and pause" parameters are ignored when this option is selected.
#' 
#' @param fssObj       Fitted model object of class "FSS".
#'                     It is the analysis results returned after calling the
#'                     \code{\link{FSS}} function.
#' @param option       Plotting option. It must be one of the following values (described above):
#'                     "full data",
#'                     "overall time series",
#'                     "specific time series",
#'                     "possible outbreaks",
#'                     "overall tree",
#'                     "sub-trees", 
#'                     "partitions", and
#'                     "outbreak detection".
#' @param cellIndices  It is only used when option = "specific time series". 
#'                     It is a mx2 matrix containing (row, col) coordinates of m cells of 
#'                     resulting grid. The time series of these m cells will be plotted.
#'                     If it is not given, the "cellIndices" specified in the 
#'                     \code{\link{FSS}} function will be used instead. 
#' @param dspSlices    A vector containing time slices to be displayed.
#' @param titles       A vector of strings containing titles of graphs for possible 
#'                     multiple figures.
#' @param pause        Whether to pause for user to press [ENTER] to continue.
#' 
#' @export
#' 
#' @seealso \code{\link{FSS}}
#' 
#' @examples
#' library(RecurPart)
#' 
#' cellIndices = rbind(c(2,9),c(2,10),c(3,9),c(3,10),
#'                     c(7,5),c(7,6),c(8,6),c(8,7))
#' 
#' dspSlices <- 1:4   # a vector containing time slices to be displayed
#' 
#' # allocate the first 24 months to find monthly means of each cell, and
#' # assuming 30 days for a month
#' test <- FSS(AEGISS_ixyt, num1Slice=30, num1Period=12, numPeriods=2, numX=10, numY=9, 
#'             plotIte=8, cellIndices=cellIndices, data_poly=AEGISS_poly, lambda=0.9, 
#'             dspSlices=dspSlices, fillThreshold=0.95, colPoints="cyan")
#' 
#' plotFSS(test, option="full data", titles=c("AEGISS all data"), pause=TRUE)
#' plotFSS(test, option="overall time series", pause=TRUE)
#' plotFSS(test, option="specific time series")
#' plotFSS(test, option="possible outbreaks")
#' plotFSS(test, option="overall tree")
#' plotFSS(test, option="sub-trees")
#' plotFSS(test, option="partitions")
#' plotFSS(test, option="outbreak detection")
#' 
plotFSS <- function(fssObj, option="full data", cellIndices=NULL, dspSlices=NULL,
										titles=NULL, pause=TRUE) {
	option = tolower(option)
	
	if (option == "full data") plotFullData(fssObj, titles)
	
	if (option == "overall time series") plotOverallTimeSeries(fssObj, titles)
	
	if (option == "specific time series") {
		plotSpecificTimeSeries(fssObj, cellIndices)
		pause = FALSE
		#return ()
	}
	
	if (option == "possible outbreaks") {
		plotPossibleOutbreaks(fssObj, dspSlices)
		pause = FALSE
		#return ()
	}
		
	if (option == "overall tree") {
		plotOverallTree(fssObj)
		pause = FALSE
		#return ()
	}
	
	if (option == "sub-trees") {
		plotSubTrees(fssObj)
		pause = FALSE
		#return ()
	}
	
	if (option == "partitions") {
		plotPartitions(fssObj)
		pause = FALSE
		#return ()
	}
	
	if (option == "outbreak detection") {
		plotOutBkDetect(fssObj)
		pause = FALSE
		#return ()
	}
	
	if (pause)
		invisible(readline(prompt="Please press [Enter] to continue"))
}

plotFullData <- function(fssObj, titles) {
	xmin = fssObj$xy$Dim$xmin; xmax = fssObj$xy$Dim$xmax
	ymin = fssObj$xy$Dim$ymin; ymax = fssObj$xy$Dim$ymax
	
	myLayout <- matrix(c(1,1), nrow=1, ncol=1, byrow=TRUE)
	nf <- layout(mat=myLayout)
	
	par(mar=c(4,2.5,2,2))
	plot(c(xmin,xmax), c(ymin,ymax), type="n", xlab="", ylab="", axes=TRUE)
	
	plotBoundaries(fssObj$poly)
	plotPoints(fssObj$xy$All,fssObj$colors$colPoints)
	
	if (is.null(titles)) title = "Showing counts of all the data"
	else title = titles[1]
	title(main=title, font.main=1, family="sans", cex.main=1.6)
	
	plotGrid(fssObj$xy$Dim, fssObj$Obs$All, fssObj$inputArgs$numY, fssObj$inputArgs$numX)
}

plotGrid <- function(Dim, Obs, numY, numX)
{
	xmin = Dim$xmin; xmax = Dim$xmax; xlen = (xmax-xmin)/numX
	ymin = Dim$ymin; ymax = Dim$ymax; ylen = (ymax-ymin)/numY
	
	# plot horizontal lines:
	x1 = xmin; x2 = xmax
	for (i in 0:numY) {
		y1 = y2 = ymin+i*ylen
		lines(c(x1,x2),c(y1,y2),col='red')
	}
	
	# plot vertical lines:
	y1 = ymin; y2 = ymax
	for (j in 0:numX) {
		x1 = x2 = xmin+j*xlen
		lines(c(x1,x2),c(y1,y2),col='red')
	}
	
	# plot text for counts:
	# each cell in (numY*numX) grid contains total counts for all data
	for (row in 1:numY)
		for (col in 1:numX)
		{
			ch <- round(Obs[row,col])
			if (ch>0) {
				y = ymin + ylen*(numY-row+1-0.5)
				x = xmin + xlen*(col-0.5)
				if (numX < 15)
					text(x,y,as.character(ch),cex=1,col='black')
				else
					text(x,y,as.character(ch),cex=0.7,col='black')
			}
		}
}

plotPoints <- function(xyPart, colPoints) {
	points(xyPart[,1], xyPart[,2], pch='.', cex=0.75, col=colPoints)
}

plotBoundaries <- function(poly) {
	n = dim(poly)[1]
	x0 = x1 = poly[1,1]
	y0 = y1 = poly[1,2]
	lineWidth = 1   # 1 or 2
	for (i in 2:n) {
		x2 = poly[i,1]; y2 = poly[i,2]
		lines(c(x1,x2),c(y1,y2),lwd=lineWidth)
		x1 = x2; y1 = y2;
	}
	lines(c(x1,x0),c(y1,y0),lwd=lineWidth)
}

