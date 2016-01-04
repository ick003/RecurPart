\name{RecurPart}
\alias{RecurPart}
\alias{RecurPart}
\docType{package}
\title{
\packageTitle{RecurPart}
}
\description{
\packageDescription{RecurPart}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{RecurPart}
\packageIndices{RecurPart}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
\packageAuthor{RecurPart}

Maintainer: \packageMaintainer{RecurPart}
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{ package }
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
