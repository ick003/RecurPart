#' @title Forward Selection Scan (FSS) Statistic
#' 
#' @description
#' Spatial-Temporal Disease Surveillance using Forward Selection Scan Statistic.
#' 
#' The approach starts by dividing the target geographical regions into a lattice.
#' Secondly it smooths the time series of lattice cell counts using multivariate 
#' exponential weighted moving averages. 
#' Thirdly, these EWMA cell counts are spatially smoothed to reduce spatial noise 
#' and leave the spatial signal. 
#' The fourth step uses forward selection approach to scanning mutually exclusive 
#' and exhaustive rectangular regions of dynamic dimensions. 
#' In the fifth step, it prunes away all insignificant scanned regions where
#' counts are not significantly higher than expected.
#' 
#' This package is under development.
#' 
#' For a complete list of exported functions, use library(help = "RecurPart")
#' 
#' @author
#' Adrien Ickowicz <\email{Adrien.Ickowicz@@csiro.au}>, 
#' Ross Sparks <\email{Ross.Sparks@@csiro.au}>,
#' Thomas Lo <\email{Thomas.Lo@@csiro.au}>
#' 
#' Maintainer: Adrien Ickowicz <\email{Adrien.Ickowicz@@csiro.au}>
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
#' @references
#' Sparks, R. S. and Ickowicz, A. (2013). Spatio-Temporal Disease Surveillance: 
#' Forward Selection Scan Statistic. arXiv preprint arXiv:1309.7721.
#' 
#' @docType package
#' @name RecurPart-package
NULL

#' @title ID numbers and space-time locations.
#'
#' @description
#' A dataset containing id numbers and space-time locations,
#' (id,x,y,t) for 10,572 cases of non-specific gastrointestinal 
#' disease in the county of Hampshire, UK, as reported to NHS Direct 
#' (a phone-in triage service operating within the UK's National 
#' Health Service). Each (x,y)-location corresponds to the centroid of 
#' the unit post-code of the residential address of the person making 
#' the call to NHS Direct. The unit of distance is 1 metre. The unit of 
#' time is 1 day, with day 1 corresponding to 1 January 2001. 
#' 
#' @format A data frame with 10572 rows and 4 variables. 
#' The variables are as follows:
#' 
#' \itemize{
#'   \item id. ID number
#'   \item x. x-coordinate of centroid
#'   \item y. y-coordinate of centroid
#'   \item t. Time in day
#' }
#'
#' @details
#' ACKNOWLEDGEMENT. Data were provided by Prof.Peter Diggle, Lancaster
#' University and Dr Peter Hawtin, Health Protection Agency,Southampton Laboratory.
#' The data were derived from anonymised data-sets collected during
#' project AEGISS, a collaborative surveillance project based in the south
#' of England.  The project was financially supported by the Food Standards Agency.
#' 
#' CITATION: Diggle, P, Knorr-Held, L, Rowlingson, B, Su, T, Hawtin, P
#' and Bryant, T (2003). On-line monitoring of public health surveillance 
#' data. In "Monitoring the Health of Populations: Statistical Principles
#' and Methods for Public Health Surveillance" editors R. Brookmeyer and
#' D.F. Stroup, pages 233-66. Oxford : Oxford University Press.
#' 
#' @source \url{http://www.lancaster.ac.uk/staff/diggle/pointpatternbook/datasets/}
#' 
#' @docType data
#' @name AEGISS_ixyt
NULL

#' @title Vertices of a 120-sided polygon.
#'
#' @description
#' A dataset containing vertices of a 120-sided polygon
#' representing the boundary of the study-region (Hampshire,UK).
#' The unit of distance is 1 kilometre.. 
#'
#' @format A data frame with 119 rows and 2 variables. 
#' The variables are as follows:
#' 
#' \itemize{
#'   \item x. x-coordinate of vertex
#'   \item y. y-coordinate of vertex
#' }
#'
#' @details
#' ACKNOWLEDGEMENT. Data were provided by Prof.Peter Diggle, Lancaster
#' University and Dr Peter Hawtin, Health Protection Agency,Southampton Laboratory.
#' The data were derived from anonymised data-sets collected during
#' project AEGISS, a collaborative surveillance project based in the south
#' of England.  The project was financially supported by the Food Standards Agency.
#' 
#' CITATION: Diggle, P, Knorr-Held, L, Rowlingson, B, Su, T, Hawtin, P
#' and Bryant, T (2003). On-line monitoring of public health surveillance 
#' data. In "Monitoring the Health of Populations: Statistical Principles
#' and Methods for Public Health Surveillance" editors R. Brookmeyer and
#' D.F. Stroup, pages 233-66. Oxford : Oxford University Press.
#' 
#' @source \url{http://www.lancaster.ac.uk/staff/diggle/pointpatternbook/datasets/}
#' 
#' @docType data
#' @name AEGISS_poly
NULL
