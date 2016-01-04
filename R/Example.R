library(RecurPart)

cellIndices = rbind(c(2,9),c(2,10),c(3,9),c(3,10),
										c(7,5),c(7,6),c(8,6),c(8,7))

dspSlices <- 1:13   # a vector containing time slices to be displayed

# allocate the first 24 months to find monthly means of each cell, and
# assuming 30 days for a month
test <- FSS(AEGISS_ixyt, num1Slice=30, num1Period=12, numPeriods=2, numX=10, numY=9, 
						plotIte=8, cellIndices=cellIndices, data_poly=AEGISS_poly, lambda=0.9, 
						dspSlices=dspSlices, fillThreshold=0.95, colPoints="cyan")

plotFSS(test, option="full data", titles=c("AEGISS all data"), pause=TRUE)

plotFSS(test, option="overall time series", titles=NULL, pause=TRUE)

plotFSS(test, option="specific time series")

plotFSS(test, option="possible outbreaks")

plotFSS(test, option="overall tree")

plotFSS(test, option="sub-trees")

plotFSS(test, option="partitions")

plotFSS(test, option="outbreak detection")
