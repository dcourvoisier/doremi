#' Meaurements of cardiac frequency in 21 patients during medical effort tests
#'
#' A dataset containing time, cardiac frequency and load of the resistive bicycle run by patients during effort tests
#'
#' @format A data frame with 1686 rows of 4 variables
#' \describe{
#'   \item{id}{arbitrary identifier of the patient, dimensionless}
#'   \item{timecol}{time since the beginning of the test, in seconds (s)}
#'   \item{load}{load of the resistive bycicle, effort that the patient needs to do, in watts (W)}
#'   \item{hr}{patient's cardiac rythm, in heart beats per minute (1/min)}
#'
#' }
#' @docType data
#' @usage data(cardio)
"cardio"
#' Meaurements of response time of 17 individuals when carrying out mental rotation tasks
#'
#' A dataset containing trial number and response time of individuals that participated in the study by Courvoisier et al.(2013)
#'(available at https://doi.org/10.1016/j.yhbeh.2012.12.007)
#'
#' @format A data frame with 619 rows of 5 variables
#' \describe{
#'   \item{id}{arbitrary identifier of the individual, dimensionless}
#'   \item{sex}{sex of the individual, as the study highlighted the difference in response time according to sex, dimensionless}
#'   \item{numtrials}{number of trial to which the mean response time belongs, dimensionless}
#'   \item{meanRT}{mean response time of the individual to execute the mental rotation task, in milliseconds (ms)}
#'   \item{logmeanRT}{natural logarithm of the mean response time, in milliseconds (ms)}
#'
#' }
#' @docType data
#' @usage data(rotation)
"rotation"
