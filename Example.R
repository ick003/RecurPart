library(RecurPart)

cellIndices = rbind(c(2,9),c(2,10),c(3,9),c(3,10),
										c(7,5),c(7,6),c(8,6),c(8,7))

#plotFlags <- c(F, F, F, F, F, F, F, T, T, T)
#plotFlags <- c(T, T, T, T, T, T, T, T, T, T)
plotFlags <- c(T, F, F, F, F, F, F, F, F, F)

dspSlices <- 1:10   # a vector containing time slices to be displayed

plotIte <- 8
# plotIte <- 5     # to test option="overall tree"

# allocate the first 24 months to find monthly means of each cell, and
# assuming 30 days for a month
test <- FSS(AEGISS_ixyt, num1Slice=30, num1Period=12, numPeriods=2, numX=10, numY=9, 
						plotIte=plotIte, cellIndices=cellIndices, data_poly=AEGISS_poly, lambda=0.9, 
						dspSlices=dspSlices, fillThreshold=0.95, colPoints="cyan")

if (plotFlags[1]) plotFSS(test, option="full data", titles=c("AEGISS all data"), pause=TRUE)

if (plotFlags[2]) plotFSS(test, option="overall time series", titles=NULL, pause=TRUE)

if (plotFlags[3]) plotFSS(test, option="specific time series")

if (plotFlags[4]) plotFSS(test, option="possible outbreaks")

if (plotFlags[5]) plotFSS(test, option="overall tree")

if (plotFlags[6]) plotFSS(test, option="sub-trees")

if (plotFlags[7]) plotFSS(test, option="partitions")

if (plotFlags[8]) plotFSS(test, option="outbreak detection")

tempCellIndices = rbind(c(3,9),c(3,10))
if (plotFlags[9]) plotFSS(test, option="specific time series", cellIndices=tempCellIndices)

tempDspSlices <- 1:3   # a vector containing time slices to be displayed
if (plotFlags[10]) plotFSS(test, option="possible outbreaks", dspSlices=tempDspSlices)
