#' Measurements of cardiac frequency in 21 patients during effort tests
#'
#' Data containing time, cardiac frequency and load of a resistive bicycle run by patients during effort tests
#'
#' @format A data frame with 1686 rows and 4 variables
#' \describe{
#'   \item{id}{positive integer, arbitrary identifier of the patient}
#'   \item{time}{positive real number, time since the beginning of the test, in seconds (s)}
#'   \item{load}{positive real number, load of the resistive bicycle, effort that the patient needs to do, in watts (W)}
#'   \item{hr}{positive real number, patient's cardiac rhythm, in heart beats per minute (1/min)}
#'
#' }
#' @source Mongin et al. 2018, under review (future DOI will be inserted here)
#' @docType data
#' @usage data(cardio)
"cardio"

#' Measurements of response time of 17 individuals when carrying out mental rotation tasks
#'
#' Data containing reaction time to a mental rotation task over a 60 day period for 17 individuals
#' \href{https://doi.org/10.1016/j.yhbeh.2012.12.007}{(Courvoisier et al., 2013)}.
#'
#' @format A data frame with 619 rows and 5 variables
#' \describe{
#'   \item{id}{positive integer, arbitrary identifier of the individual}
#'   \item{sex}{character, sex of the individual, as the study highlighted the difference in response time according to sex}
#'   \item{days}{positive integer,day since the beginning of the experiment}
#'   \item{meanRT}{positive integer, mean response time of the individual to execute the mental rotation task, in milliseconds (ms)}
#'   \item{logmeanRT}{natural logarithm of the mean response time}
#'
#' }
#' @source Delphine S. Courvoisier, Olivier Renaud, Christian Geiser, Kerstin Paschke, Kevin Gaudy, Kirsten Jordan,
#' Sex hormones and mental rotation: An intensive longitudinal investigation,
#'
#' Hormones and Behavior,
#'
#' Volume 63, Issue 2,
#'
#' 2013,
#'
#' Pages 345-351,
#'
#' https://doi.org/10.1016/j.yhbeh.2012.12.007

#' @docType data
#' @usage data(rotation)
"rotation"
