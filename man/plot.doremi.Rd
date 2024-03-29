% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/others.R
\name{plot.doremi}
\alias{plot.doremi}
\title{S3 method to plot DOREMI objects}
\usage{
\method{plot}{doremi}(x, ..., id = NULL)
}
\arguments{
\item{x, }{DOREMI object resulting from \code{\link{analyze.1order}} or \code{\link{analyze.2order}} analysis}

\item{...}{includes the additional arguments inherited from the generic plot method}

\item{id}{Identifiers of the individuals to be represented in the plot.
By default, it will print the first six individuals.}
}
\value{
Returns a plot with axis labels, legend and title. The axis labels and legend include the names of the variables set as input arguments.
The title includes the name of the DOREMI object result of the analysis. The function uses \code{\link[ggplot2]{ggplot}}
to generate the graphs and so it is possible to override the values of axis labels, legend and title through ggplot commands.
}
\description{
\code{plot.doremi} generates a plot with the observed values of the signal, the excitation values and the fitted
signal over time for each individual.
}
\examples{
mydata <- generate.panel.1order(time= 0:100,
                                excitation = sin(0:100),
                                y0 = 0,
                                t0 = 0,
                                tau = 2,
                                k = 1,
                                yeq = 0,
                                nind = 2,
                                internoise = 0.1,
                                intranoise = 8)
myresult <- analyze.1order(data = mydata,
                           id = "id",
                           input = "excitation",
                           time = "time",
                           signal = "signal")
plot(myresult)
}
