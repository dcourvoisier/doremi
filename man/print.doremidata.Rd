% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/others.R
\name{print.doremidata}
\alias{print.doremidata}
\title{S3 method to print DOREMI data objects}
\usage{
\method{print}{doremidata}(x, ...)
}
\arguments{
\item{x}{DOREMIDATA object}

\item{...}{includes the additional arguments inherited from the generic print method}
}
\value{
Returns the table $data of the DOREMIDATA object
}
\description{
\code{print.doremidata} prints the most important results of a  DOREMIDATA object
}
\examples{
mydata <- generate.panel.remi(nind = 5,
                           dampingtime = 10,
                           amplitude = c(5,10),
                           nexc = 2,
                           duration = 20,
                           deltatf = 2,
                           tmax = 200,
                           minspacing = 0,
                           internoise = 0.2,
                           intranoise = 0.1)
mydata
}
