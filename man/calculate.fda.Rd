% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doremi.R
\name{calculate.fda}
\alias{calculate.fda}
\title{Calculation of derivatives using the Functional Data Analysis (FDA) method.}
\usage{
calculate.fda(signal, time, spar = NULL, order = NULL)
}
\arguments{
\item{signal}{is a vector containing the data from which the derivative is estimated.}

\item{time}{is a vector containing the time values corresponding to the signal. Arguments signal and time must have the same length.}

\item{spar}{is the smoothing parameter used by the roughness penalty function in the smooth.spline R function.}

\item{order}{parameter not used, for consistency with calculate.glla and calculate.gold}
}
\value{
Returns a list containing two elements:

dtime- contains the initial time values provided.

dsignal- is a data.frame containing three columns and the same number of rows as the signal.
The first column is the signal data points, the second is the first derivative evaluated at those points,
and the third is the second derivative evaluated at those points.
}
\description{
\code{calculate.fda} estimates the derivatives of a variable using the FDA
method described in several sources, such as in \doi{10.1007/b98888}{Ramsay et al. (2009)}
and  \doi{10.1080/00273171.2015.1123138}{Chow et al. (2016)}.
This method estimates a spline function that fits all the data points and then derivates this function to estimate derivatives at those points.
In order for the derivatives to exist, the function must be smooth. A roughness penalty function controlled by a smoothing parameter is then used.
The estimations are done by using the R's base smooth.spline function.
}
\examples{
#In the following example the derivatives for the function y(t) = t^2 are calculated.
#The expected results are:
#y'(t) = 2t and y''(t) = 2
time <- c(1:500)/100
signal <- time^2
result <- calculate.fda(signal = signal, time = time)

}
\keyword{derivative}
\keyword{fda}
\keyword{spline}
