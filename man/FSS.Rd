\name{FSS}
\alias{FSS}

\title{Forward Scan Statistic}
\description{Perform a forward scan statistic analysis.}

\usage{
FSS(...)

\method{FSS}{default}(data_xyt,...)
}

\arguments{
\item{data_xyt}{ a spatial-temporal point process data frame.}
\item{num1Slice}{ number of time points in one slice.}
\item{num1Period}{ number of slices in one period.}
\item{numX}{ number of cells on the x-axis.}
\item{numY}{ number of cells on the y-axis.}
\item{tempSmooth}{ logical. Whether or not apply a temporal smoothing to the data.}
\item{spatSmooth}{ logical. Whether or not apply a spatial smoothing to the data.}
\item{alpha}{ value of the temporal smoother.}
\item{lambda}{ value of the spatial smoother.}
\item{num1Day}{ number of time points per day.}
\item{\dots}{ Additional parameters.}
}

\value{
An object of class list, including elements
\item{xy}{ original point process.}
\item{poly}{ polygon of the spatial area (if specified).}
\item{Obs}{ transformed data to feed in the algorithm.}
\item{cells}{ cells of interest (if specified).}
\item{inputArgs}{ list of the input arguments.}
\item{originalSlices}{ cell-repartition of the data.}
\item{smoothedSlices}{ smoothed cell repartition of the data.}
\item{partResults}{ dataframe of all the iterations of the algorithm, with the partitioning decision and the offsprings.}
\item{colors}{ set of chosen colors for plotting.}

}

\author{
  \packageAuthor{RecurPart}
  
  Maintainer: \packageMaintainer{RecurPart}
}

\examples{
cellIndices = rbind(c(2,9),c(2,10),c(3,9),c(3,10),
                    c(7,5),c(7,6),c(8,6),c(8,7))

dspSlices <- 1:13   # a vector containing time slices to be displayed

# allocate the first 24 months to find monthly means of each cell, and
# assuming 30 days for a month

test.fss <- FSS(AEGISS_ixyt, num1Slice=30, num1Period=12, 
            numPeriods=2, numX=10, numY=9, plotIte=8, 
            cellIndices=cellIndices, data_poly=AEGISS_poly,
            lambda=0.9, dspSlices=dspSlices, fillThreshold=0.95,
            colPoints="cyan")

}

\keyword{FSS}


