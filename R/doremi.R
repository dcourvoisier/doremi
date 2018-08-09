# calculate.gold ----------------------------------------------------------
#' Calculation of derivatives using the GOLD method
#'
#' \code{calculate.gold} estimates the derivatives of a variable using the Generalized Orthogonal Local Derivative (GOLD)
#' method described in \href{https://doi.org/10.1080/00273171.2010.498294}{Deboeck (2010)}.
#' This method allows calculating over a number of measurement points (called the embedding number) the first derivative with errors uncorrelated with the signal.
#' It was generalized for non-equidistant time points (variable time steps), in order to account for missing observations.
#' @param signal is a vector containing the data from which the derivative is estimated.
#' @param time is a vector containing the time values corresponding to the signal. Arguments signal and time must have the same length.
#' @param embedding is an integer indicating the number of points to consider for derivative calculation. Embedding must be greater than 1 because at least
#' two points are needed for the calculation of the first derivative and at least 3 for the calculation of the second derivative.
#' @keywords derivative, embed, rollmean
#' @return Returns a list containing three columns:
#'
#' dtime- contains the time values in which the derivative was calculated. That is, the moving average of the input time over embedding points.
#'
#' dsignal- is a data.frame containing three columns and the same number of rows as the signal.
#' The first column is the moving average of the signal over embedding points, the second is the first derivative,
#' and the third is the second derivative.
#'
#' embedding- contains the number of points used for the derivative calculation, which is constant.
#'
#' @examples
#' #In the following example the derivatives for the function y(t) = t^2 are calculated.
#' #The expected results are:
#' #y'(t) = 2t and y''(t) = 2
#' time <- c(1:500)/100
#' signal <- time^2
#' result <- calculate.gold(signal = signal, time = time, embedding = 5)
#'
#'@export
#'@importFrom zoo rollmean
calculate.gold <-  function(signal,
                            time,
                            embedding = 2){
  #Error management
  if (length(signal) != length(time)){
    stop("signal and time vectors should have the same length.\n")
  }
  if (length(signal) <= embedding){
    stop("Signal and time vectors should have a length greater than embedding.\n")
  }
  if (embedding > 2){
    #tembed (time embedded) is a matrix, containing in each line the groups of
    #"n" points with n being the embedding value from which to calculate the
    #roll means For instance if time=1:10 and embedding=3, the matrix tembed
    #will have: row1: 1 2 3; row2: 2 3 4...row8:8 9 10).
    tembed <- embed(time, embedding)
    #The "stats" library "embed" function provides a matrix in which the groups
    #of values are inversed by column and thus the next operation is to format
    #these values so that the matrix contain the groups by line and from left to
    #right.
    tembed <- tembed[, ncol(tembed):1]

    #Creation of the D matrix for estimation of derivatives up to second
    #derivative. According to the paper mentioned, equations (11) and (12) this
    #is "...the diagonal matrix with scaling constants that convert the
    #polynomial estimates to derivative estimates".
    D <- cbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 0.5))
    #Cbind does a column binding of the vectors


    #Creation of the empty derivative matrix Lines: It will have as many lines
    #as groups of time points the embed function gives (thus, the number of lines
    #of the tembed matrix) Columns: It will have 3 columns because: First column
    #will contain the signal mean value in the time points considered in
    #"embedding". Second column will contain the signal derivative in those same
    #time points. Third column will contain the signal second derivative in those
    #same time points
    derivative <- matrix(NA, nrow = nrow(tembed), ncol = 3)


    #Xembed is a matrix containing the values of the signal in the time points
    #considered in "embedding" (no mean calculated yet) As for tembed, it
    #contains in each line the groups of "n" points of the function with n being
    #the embedding value to PREPARE for the calculation of the roll mean
    Xembed <- embed(signal, embedding)
    Xembed <- Xembed[, ncol(Xembed):1]

    for (k in 1:nrow(tembed)){
      # Loop repeated in each tembed line (each group of time values)
      # To take into account variable delta t.
      t <- tembed[k, ] - tembed[k, (embedding + 1) / 2]
      #time vector, containing deltat, centered in 0
      E <- matrix(NA, nrow = 3, ncol = length(t))
      #Construction of the empty Theta matrix, equation (10) But will be filled
      #from the polynomials calculated in equation (9). It has 3 rows as we are
      #calculating up to the second derivative. The two for loops resolve the
      #paper's equation (9) (polynomial's coefficients)
      for (i in 1:length(t)){
        E[1, i] <- 1
        E[2, i] <- t[i]
      }
      for (i in 1:length(t)){
        E[3, i] <- t[i] ^ 2 - sum(t ^ 2) / (sum(E[1, ])) - t[i] * sum(t ^ 3) / sum(t ^ 2)
      }
      # And with E (Theta) and D it is possible to calculate the orthogonal
      # matrix W, equation (14)
      L <- D %*% E
      W <- t(L) %*% solve(L %*% t(L))
      #And with this matrix, solve the differential equation as Y=X*W
      derivative[k, ] <- Xembed[k, ] %*% W
      #Derivative calculated in the time row k
    }
    derivative <- rbind(derivative, matrix(data = NA, ncol = 3, nrow = embedding - 1))
    # Addition of NA so that the derivative rows are the same as those of signal

  } else if (embedding == 2){
    #warning("Only first derivative can be calculated with an embedding of 2.\n")
    derivative <- cbind(rollmean(signal, embedding), diff(signal) / diff(time))
    #Appends the two columns: 1. mean time values and 2. The span calculated as
    #the difference of two signal values (going forward) divided by the time
    #interval
    derivative <- rbind(derivative,matrix(data = NA, ncol = 2, nrow = embedding - 1))
    # Addition of NA so that the derivative rows are the same as those of signal
  } else{
    stop("Embedding should be >=2 for the calculation of derivatives.\n")
  }
  time_derivative <-c(rollmean(time, embedding), rep(NA, embedding - 1))
  # Addition of NA so that the time rows are the same as those of signal

  returnobject <- list("dtime" = time_derivative,
                       "dsignal" = derivative,
                       "embedding" = embedding)
  return(returnobject)
}
# generate.excitation -----------------------------------------------------
#' Excitation signal generation
#'
#' \code{generate.excitation} generates a vector of randomly located square pulses
#' with a given amplitude, duration and spacing between the pulses. A pulse is where the excitation passes from value 0 to value amplitude
#' for a given duration and then returns back to 0, thus producing a square shape.
#' @param amplitude  is vector of values different than 0 indicating the amplitude of the excitation. It should contain as many values
#' as the number of pulses (nexc). If the elements are less than the number of pulses, the amplitude vector will be "recycled" and the elements from it will be repeated until
#' all the pulses are covered (for instance, if the number of excitations nexc is 6 and the amplitude vector has two elements, pulses 1,3 and 5 will
#' have the same amplitude as the first element of the amplitude vector and pulses 2,4 and 6 that of the second element).
#' @param nexc is an integer greater than 0 indicating the number of pulses to generate.
#' @param duration is a vector of values greater than 0 indicating the duration of each pulse in time units. It should have as many elements as the number of pulses (nexc). If
#' the elements are less than the number of pulses, the amplitude vector will be "recycled" and the elements from it will be repeated until
#' all the pulses are covered.
#' @param deltatf is a value greater than 0 indicating the time step between two consecutive data points.
#' @param tmax is a value greater than 0 indicating the maximum time range of the excitation vector in time units. The time vector generated will go from 0 to tmax.
#' @param minspacing as pulses are generated randomly, minspacing is a value greater than 0 that indicates minimum spacing between pulses, thus avoiding
#' overlapping of the pulses in time. It can be 0 indicating that pulses can follow one another.
#' @keywords excitation, simulation
#' @return Returns two vectors:
#'
#' E- vector containing the values of the excitation generated.
#'
#' t- vector containing the values of time generated.
#' @details Used for simulations in the context of the package. Beware that the following condition should apply:
#' \deqn{tmax >= (duration+minspacing)*nexc}
#' so that the pulses "fit" in the time lapse defined.
#' Compared to pulsew from the seewave package this function can generate pulses of different duration and amplitude.
#' @examples
#' generate.excitation (amplitude = 3,
#'                      nexc = 6,
#'                      duration = 2,
#'                      deltatf = 1,
#'                      tmax = 200,
#'                      minspacing = 2)
#' #Vector of length 201 (deltatf x tmax + 1 as it includes 0 as initial time value)
#' generate.excitation (amplitude = c(1,10,20),
#'                      nexc = 3,
#'                      duration = c(1,2,4),
#'                      deltatf = 0.5,
#'                      tmax = 100,
#'                      minspacing = 10)
#'@export
#'@importFrom utils head
# Excitation signal generation.
generate.excitation = function(amplitude = 1,
                               nexc = 1,
                               duration = 2,
                               deltatf = 0.1,
                               tmax = 10,
                               minspacing = 1)
{
  #Error management
  if (any(duration <= 0)){
    stop("Invalid input parameters. Pulse duration must be greater than 0.\n")
  }
  if (any(amplitude == 0)){
    stop("Invalid input parameters. Amplitude must be different from 0.\n")
  }
  if (nexc <= 0){
    stop("Invalid input parameters. At least one excitation must be defined.\n")
  }
  if (nexc < length(duration) | nexc < length(amplitude)){
    stop("The number of excitations nexc is smaller than the number of elements in amplitude and/or duration.\n")
  }

  if (nexc > length(duration)){
    if (length(duration) > 1){
      warning("The number of excitations nexc was higher than the durations defined. Durations were recycled.\n")
    }
    duration <- rep(duration, ceiling(nexc/length(duration)))
    duration <- duration[1:nexc]
  }
  if (nexc > length(amplitude)){
    if (length(amplitude) > 1){
      warning("The number of excitations nexc was higher than the amplitudes defined. Amplitudes were recycled.\n")
    }
    amplitude <- rep(amplitude,ceiling(nexc/length(amplitude)));
    amplitude <- amplitude[1:nexc]
  }
  if(tmax < (sum(duration) + minspacing * (nexc-1))){stop("Invalid input parameters. tmax should be greater than (duration + minspacing) * nexc.\n")}

  #Generation of time vector
  tim <- seq(0, tmax, deltatf)

  found <- FALSE #indicates final distribution of pulses has not been found yet
  while (!found){
    tal <- c(sort(sample(0:tmax, nexc, replace = F)), tmax) #initial sampling. tal is a vector of time values in which the pulses will start
    if(all(diff(tal) >= (minspacing + duration))){found <- TRUE} #means all the pulses fitted and we exit the loop
  }

  #Generation of excitation
  E <- rep(0,length(tim)) #initialize excitation vector with 0

  for (i in 1:nexc){
    E[(tim >= tal[i]) & (tim <= (tal[i] + duration[i]))] <- amplitude[i] #fill vector with the amplitudes for the corresponding pulse in the lapses
    #of time defined by tal and tal+duration
  }

  data <- list(exc = E, t = tim)
  return(data)
}
# generate.remi ----------------------------------------------------
#' Generation of the first order differential equation solution
#'
#' \code{generate.remi} returns a vector containing the convolution of the Green function
#' of the first order differential equation with a given damping time and an excitation term.
#'
#' @param dampingtime Signal damping time. It represents the characteristic response time of the solution of the differential equation.
#' It should be positive.
#' @param inputvec Is a vector containing the values of the excitation signal.
#' @param inputtim Is a vector containing the time values corresponding to the excitation signal.
#' @keywords Green, excitation, first order differential equation
#' @return Returns a list containing two elements:
#' \itemize{
#'   \item  y is a vector containing the values calculated from the convolution of the Green function and the excitation vector.
#'   \item  t is a vector containing the corresponding time values
#' }
#' @examples
#' generate.remi(10,rep(c(0,1,0),c(20,30,50)),1:100)
#'@export
#'@importFrom stats convolve
generate.remi <- function(dampingtime,
                          inputvec,
                          inputtim){
  #Error management
  #If inputvec is a scalar, the function warns the user that it should be a vector containing the values of the excitation signal
  if (length(inputvec) <= 1 | length(inputvec) != length(inputtim)) {
    stop("Both the excitation (inputvec) and its time values (inputtim) should be vectors and have the same length.\n")
  }
  #if inputvec is a character or a matrix, the function stops
  if (is.matrix(inputvec) | is.character(inputvec)) stop("Inputvec should be a vector.")

  green <- exp(-1 / dampingtime * (inputtim - inputtim[1])) #Calculation of green function values.inputtim vector is forced to start in 0

  #as convolution doesn't take time into account, only number of points, a coefficient is needed to reduce amplitude
  convol <- (max(inputtim) - min(inputtim)) / (length(inputtim) - 1) * rev(convolve(rev(green), inputvec, type = "open"))
  convol <- convol[1:length(inputvec)]

  return(list(y = convol, t = inputtim))              #Returns the result of the convolution in the points corresponding to the excitation time vector.
}
# remi --------------------------------------------------------------------
#' DOREMI first order analysis function
#'
#' \code{remi}  estimates the coefficients of a first order differential equation of the form:
#' \deqn{\frac{1}{\gamma} \dot{y}(t) = - y(t) + \epsilon E(t) + eqvalue}
#' using linear mixed-effect models.
#' Where y(t) is the individual's signal, \eqn{\dot{y}(t)} is the derivative and E(t) is the excitation.
#' @param data Is a data frame containing at least one column, that is the signal to be analyzed.
#' @param id Is a STRING containing the NAME of the column of data containing the identifier of the individual.
#' If this parameter is not entered when calling the function, a single individual is assumed and a linear regression is done instead
#' of the linear mixed-effects regression.
#' @param input Is a STRING or a VECTOR OF STRINGS containing the NAME(s) of data column(s) containing the EXCITATION vector(s).
#' If this parameter is not entered when calling the function,
#' the excitation is assumed to be unknown. In this case, the linear mixed-effect regression is still carried out but no coefficient is calculated
#' for the excitation term. The function then uses the parameters estimated by the regression to carry out an exponential fit of the signal
#' and build the estimated curve.
#' The function will consider an excitation variable each column of data whose name is contained in the input vector.
#' The function will return a coefficient for each one of the excitation variables included.
#' @param time Is a STRING containing the NAME of the column of data containing the time vector. If this parameter is not entered when calling the function,
#' it is assumed that time steps are of 1 unit and the time vector is generated internally in the function.
#' @param signal Is a STRING containing the NAME of the column of the data frame containing the SIGNAL to be studied.
#' @param embedding Is a positive integer containing the number of points to be used for the calculation of the derivatives. Its value by default is 2 as at
#' least two points are needed for the calculation of the first derivative.
#' @keywords analysis, first order, exponential
#' @return Returns a summary of the fixed components for the three coefficients: damping time, excitation coefficient and equilibrium value.
#' @details The analysis performs the following linear mixed-effects regression:
#' \deqn{y_{ij}'  \sim   b_{0} +b_{0j}+b_{1} y_{ij}+b_{2} E_{ij}+u_{1j} y_{ij}+u_{2j} E_{ij}+e_{ij}}
#' with i accounting for the time and j for the different individuals. \eqn{e_{ij}} are the residuals,
#' \eqn{y_{ij}'} is the derivative calculated on embedding points and
#' y and E are the signal and the excitation averaged on embedding points.
#' The coefficients estimated to characterize the signal are calculated as follows:
#' \itemize{
#'   \item Damping time:  \eqn{\tau _{j} =  \frac{1}{ \gamma _{j} }}  with \eqn{\gamma _{j} =  b_{1} + u_{1j} }
#'   \item Excitation coefficient: \eqn{\epsilon = \frac{b_{2} + u_{2j}}{\gamma _{j}}}. It is the proportionality between the excitation and the
#'   difference between the maximum value reached by the signal and its initial value.
#'   \item Equilibrium value: \eqn{eqvalue = \frac{b_{0} + b_{0j}}{\gamma _{j}}}. It is the stable value reached in the absence of excitation.
#' }
#' The estimation is performed using the function lmer if there are several individuals or lm if there is only one.
#' With the above estimated parameters, the estimated signal can be reconstructed for
#' each individual by first performing the convolution of the excitation with the Green function using the estimated damping rate
#' and then offsetting the resulting signal with the equilibrium value.
#' The function returns five objects:
#' \enumerate{
#'  \item data- A data.frame including the input data and intermediate calculations used to prepare the variables for
#'   the fit:
#'
#'     signal_rollmean - calculation of the moving average of the signal over embedding points.
#'
#'     signal_derivate1 - calculation of the first derivative of the signal with the GOLD method in embedding points.
#'
#'     time_derivate - calculation of the moving average of the time vector over embedding points.
#'
#'     input_rollmean - calculation of the moving average of the excitation vector over embedding points.
#'  \item resultid- A data.frame including for each individual, listed by id number, the damping time, the excitation coefficient and the
#'  equilibrium value (see variables presented in the Details section).
#'  \item resultmean- A data.frame including the fixed effects of the three coefficients mentioned above.
#'  \item regression- A list containing the summary of the linear mixed-effects regression.
#'  \item estimated- A data.frame containing the estimated signal calculated as the convolution of the Green function with the
#'  estimated damping time, excitation coefficient and equilibrium value and the excitation
#'  vector with an added offset (see above). There are two extreme cases in the generation of the signal and these depend on sampling. The excitation
#'  vector is expanded to generate a pseudo-continuous signal and increase accuracy when calculating the convolution. Missing data in the exitation signal
#'  is completed by using the previous known value (exc_min in the data.frame) or using the next known value
#'  (exc_max in the data.frame). With these two inputed excitation vectors, the two extreme cases of the
#'  mean estimated signal are calculated by carrying out the convolution between the Green function (decreasing exponential with the damping time calculated by
#'  the linear mixed-effects regression). These excitation variables are then used to generate the ymin and ymax signals respectively.
#'
#'  As seen in the Description section, the print method by default prints only the resultmean element. Each one of the other objects
#'  can be accessed by indicating $ and their name after the result, for instance, for a DOREMI object called "result", it is possible
#'  to access the regression summary by typing result$regression.
#'  \item embedding - containts the embedding number used to generate the results (same as function input argument)
#' }
#' @seealso \code{\link{calculate.gold}} to compute the derivatives, for details on embedding.
#' @examples
#' myresult <- remi(data = cardio,
#'                  id = "id",
#'                  input = "load",
#'                  time = "time",
#'                  signal = "hr",
#'                  embedding = 5)
#'@export
#'@import data.table
#'@importFrom lme4 lmer
#'@importFrom lme4 lmerControl
#'@importFrom lme4 ranef
#'@importFrom stats embed
#'@importFrom stats lm
#'@importFrom zoo rollmean
#'@importFrom zoo na.locf
remi <- function(data,
                 id = NULL,
                 input = NULL,
                 time = NULL,
                 signal,
                 embedding = 2){
  #If it is a doremidata object, it takes by default the $data table from it.
  #Otherwise it takes whaterever is assigned to the data variable.
  if (first(class(data)) == "doremidata"){intdata <- copy(data$data)
  }else{intdata <- copy(data)} # Makes a copy of data so that it can rename columns freely if needed

  intdata <- setDT(intdata) # Converts data to data.table.
  noinput <- FALSE #Flag that will allow to differentiate if there is an excitation term or not when doing the regression

  intdata <- intdata[,.SD,.SDcols = c(id, input, time, signal)] # Extracting only relevant columns
  # Error management --------------------------------------------------------
  if (!is.null(id)){errorcheck(intdata, id)}
  if (!is.null(input)){errorcheck(intdata, input)
  }else{
    #If no input argument is provided, a warning will be generated indicating that the input has been set to 0.
    intdata[, input := 0]
    input = "input"
    noinput <- TRUE #This flag will be needed later as if there is no input, coefficients for the excitation term will not be calculated in the regression
    warning("No excitation signal introduced as input. Input was set to 0.\n")
  }
  if (!is.null(time)){errorcheck(intdata, time)
  }else{
    #If no time is provided, a warning will be generated and a time column with 1 time unit increment will be generated.
    intdata[, time := 1L * c(1:.N), by = id]
    time <- "time" # if no time set it to a 1 sec step vector
    warning("No time vector introduced as input. A 1 unit increment time vector was generated.\n")
  }

  #Verifying column names repeated in data table.
  if(any(duplicated(colnames(intdata)))){stop("Input datatable contains duplicated column names.\n")}

  #Suppress NA elements if there is an NA in time or in signal
  #This suppresses the entire row of the intdata table
  intdata <- intdata[!is.na(time) & !is.na(signal)]

  #If in the intdata left, there are some NA values in the excitation vector, set them to zero.
  #This is to avoid losing intdata from signal and time vectors, as sometimes when people fill in the signal intdata they put NA to mean no excitation, thus 0.
  na.to.0 <- function(x){x[is.na(x)] <- 0; x}
  intdata[, (input) := lapply(.SD, na.to.0), .SDcols = input]

  # Regression for SEVERAL individuals --------------------------------------


  # First order derivative equation mixed regression
  if (!is.null(id)){ #If the id column is not null then it is assumed that there are SEVERAL INDIVIDUALS
    print("Status: panel data")
    #Saving original names and then renaming data table to internal names
    originalnames <- c(id, time, signal, input)
    doremiexc <- paste0("input", seq(input)) # Doremi excitation vector ("input1","input2","input3"...)
    doreminames <- c("id_tmp", "time", "signal", doremiexc)
    setnames(intdata, originalnames, doreminames)

    #Rename id column to "id_tmp" and create an extra column "id" that is numeric
    #Regression works better with numeric id (comes from definition of lmer)
    intdata[, id := rep(seq(unique(id_tmp)), intdata[, .N, by = id_tmp]$N)]

    #Calculation of the signal rollmean and first derivative of the signal column
    #Paste is only used to generate new column names based on the orignal ones (concatenate strings)
    intdata[, signal_rollmean := calculate.gold(signal, time, embedding)$dsignal[, 1], by = id]
    intdata[, signal_derivate1 := calculate.gold(signal, time, embedding)$dsignal[, 2], by = id]
    intdata[, time_derivate := calculate.gold(signal, time, embedding)$dtime, by = id]

    #Calculation of the roll mean of the excitation columns if there is at least one input column
    if (!noinput){
      myfun <- function(x){x[] <- c(rollmean(x, (embedding)), rep(NA, embedding - 1)); x}
      intdata[, (paste0(doremiexc,"_rollmean")) := lapply(.SD, myfun), .SDcols = doremiexc, by = id]
    }

    #Linear mixed-effect regression
    if (noinput){ # if there is no excitation signal
      model <- tryCatch({lmer(signal_derivate1 ~ signal_rollmean + (1 + signal_rollmean | id),
                              data = intdata, REML = TRUE, control = lmerControl(calc.derivs = FALSE, optimizer = "nloptwrap"))}, error = function(e) e)
      print("Status: Unknown excitation. Linear mixed-effect model calculated.")
    }else{ # if there is one OR SEVERAL excitation signals
      model <- tryCatch({lmer(paste0("signal_derivate1 ~ signal_rollmean + (1 + ", paste(doremiexc, "rollmean ", collapse = "+", sep = "_"),
                              " + signal_rollmean | id) + ", paste(doremiexc, "rollmean ", collapse = "+",sep = "_")),
                              data = intdata, REML = TRUE, control = lmerControl(calc.derivs = FALSE, optimizer = "nloptwrap"))}, error = function(e) e)
      print("Status: One or several excitations. Linear mixed-effect model calculated.")
    }
    if (!inherits(model,"error")){ # if the regression worked
      print("Status: Linear mixed-effect model had no errors.")
      summary <- summary(model) # Summary of the regression
      random <- ranef(model)$id # Variation of the estimate coefficient over the individuals
      regression <- list(summary, random) # list to output both results: summary, and the table from ranef

      #Create tables with results
      #The first table contains the input data

      #The second table contains the mean values for gamma and thao for all the individuals (single line)
      #Generate mean results with convergence criterions
      resultmean <- setDT(list(id = "All"))

      # calculate the damping time for all signal columns: -1/damping_coeff
      resultmean[, dampingtime := -1L/summary$coefficients["signal_rollmean", "Estimate"]]

      # Extract the intercept coeff (equilibrium value)
      resultmean[, eqvalue := summary$coefficients["(Intercept)","Estimate"] * resultmean[, dampingtime]]

      #The third table contains the results for gamma and thao for each individual (one line per individual)
      resultid <- setDT(list(id = unique(intdata$id), id_tmp = unique(intdata$id_tmp)))
      setkey(resultid, id) #sorts the data table by id

      # Damping time in resultid ------------------------------------------------
      resultid[, dampingtime := resultmean[, dampingtime] + random[.GRP, "signal_rollmean"], by = id]

      # Extract the intercept (equilibrium value) calculated for each individual (present in random, regression table)
      # Offset in resultid ------------------------------------------------------
      resultid[, eqvalue := (summary$coefficients["(Intercept)", "Estimate"] + random[.GRP, "(Intercept)"]) * resultid[.GRP, dampingtime], by = id]

      # Fifth table (after regression) contains the estimated signal
      # Generation of the estimated signal for all id using remi generate FOR SEVERAL INDIVIDUALS
      # Write a warning if any of the damping times calculated was negative
      if (any(is.na(resultid[, dampingtime])) | any(resultid[, dampingtime] < 0)){
        warning("Some of the damping times calculated were negative and thus, the estimated signal was not generated for these.
                Damping times can be negative for some individuals for the following reasons: 1. The signal of
                the individual doesn't go back to equilibrium. 2.The linear mixed-effects model estimating the random
                effect showed some error messages/warnings. 3.Model misspecification.\n")
      }
      if (noinput){ #There is no excitation signal as input: exponential fit
        print("Status: Unknown excitation. Calculation of estimated signal through exponential fit.")
        #Exponential model is assumed when no input is provided and thus a new fit is necessary BY INDIVIDUAL
        #log(y-B) = gamma*t + log(A) --> LINEAR EQUATION
        #B is known, it is the intercept*thao, from the model fit with the derivative inside the analysis function
        #A is unknown. It can be found by fitting log(y-B)~ t
        #We'll call the coefficients resulting from this new fit Ap and Bp: Ap=gamma, Bp=log(A)
        intdata[, expfit_A := lm(log(abs(signal - resultid[.GRP, eqvalue])) ~ time)$coefficients[1], by = id]
        intdata[, signal_estimated :=
                  if(!is.na(resultid[.GRP, dampingtime]) && resultid[.GRP, dampingtime] > 0){
                    # if there is a damping time that has been calculated and if it is greater than 0 (decreasing exponential)
                    #Then it is assumed that the signal follows a decreasing exponential: y = A* exp(gamma*t)+B
                    #expmodel comes from fitting log(y-B)~t. Calculated above
                    exp(expfit_A) * exp(-1L / resultid[.GRP, dampingtime] * time) + resultid[.GRP, eqvalue]
                  }else{NaN}, by = id]

        #Adding calculated columns to table "estimated" and removing them from table "intdata"
        estimated <- intdata[, c("id_tmp", "time", "signal_estimated"), with = FALSE]
        intdata <- intdata[, c("expfit_A", "signal_estimated") := NULL]

      }else{
        print("Status: One or several excitation terms. Calculation of estimated signal through convolution.")
        # Extract the excitation coeff for each excitation
        intdata[, totalexc := 0]
        for (i in 1:length(input)){ #For loop to go through all the inputs
          resultmean[, paste0(doremiexc[i],"_coeff") := summary$coefficients[paste0(doremiexc[i], "_rollmean"), "Estimate"] * resultmean[, dampingtime]]

          #If variation of the excitation coefficient across individuals needed:
          #And for each individual: the mean coeff (sumary$coeff) + the variation per Individual (in random)
          #Excitation coefficient in resultid --------------------------------------
          resultid[, paste0(doremiexc[i],"_coeff") := (summary$coefficients[paste0(doremiexc[i], "_rollmean"), "Estimate"] +
                                                         random[.GRP,paste0(doremiexc[i], "_rollmean")]) * resultid[.GRP, dampingtime], by = id]
          intdata[, totalexc := totalexc + (summary$coefficients[paste0(doremiexc[i], "_rollmean"), "Estimate"] +
                                              random[.GRP,paste0(doremiexc[i], "_rollmean")]) * get(doremiexc[i]), by = id]
        }
        #The estimated table contains expanded time vector, minimum and maximum generated signals (for the two extreme scenarios of expanded excitation)
        #Generation of the expanded excitation vector according to desired deltat. Minimum and maximum values
        #Time vector and excitation vectors min and max
        intdata[, deltat := 0.1 * min(diff(time)), by = id]

        #Expanded time vector: takes the minimum and the maximum of all the time intervals and creates the vector with deltat time intervals
        estimated <- intdata[, list(time = seq(floor(min(time, na.rm = T)), ceiling(max(time, na.rm = T)), deltat[1])), by = c("id_tmp", "id")]

        #Finding the values of the original time vector in the expanded time vector by using the function "findInterval"
        IDvec <- unique(intdata[[id]])

        for (idx in seq(IDvec)){ #It is necessary to do a loop because findInterval finds the index in which the value is found in the original time vector
          #And it will be necessary to shift these indexes according to in which id we are
          leftidx <- findInterval(intdata[id == IDvec[idx], time], estimated[id == IDvec[idx], time])
          #tmpsignal contains the value of exctotal were the original times are found in the new time vector
          #and NA in the rest of the positions
          estimated[nrow(estimated[id < IDvec[idx]]) + leftidx, tmpsignal := intdata[id == IDvec[idx], totalexc]]
        }

        #The na.locf function (last observation called forward) will repeat the last non NA value.
        #We use it to repeat the values in the excitation function to the left and to the right
        estimated[, exc_min := na.locf(tmpsignal, na.rm = F, fromLast = F), by = id]
        estimated[, exc_max := na.locf(tmpsignal, na.rm = F, fromLast = T), by = id]
        estimated <- estimated[!is.na(exc_min) & !is.na(exc_max)]

        #Calculating the convolution
        estimated[, ymin := generate.remi(resultid[.GRP, dampingtime], exc_min, time)$y + resultid[.GRP, eqvalue], by = id]
        estimated[, ymax := generate.remi(resultid[.GRP, dampingtime], exc_max, time)$y + resultid[.GRP, eqvalue], by = id]

        #Removing id and tmpsignal from estimated and totalexc and total_exc from intdata
        estimated[, c("id", "tmpsignal") := NULL]
        intdata[, c("deltat","totalexc") := NULL]

      }
      #Renaming columns in $resultid, $resultmean and $estimated objects to original names

    }else{ # if the regression didn't work, a warning will be generated and tables will be set to NULL
      print("Status: Linear mixed-effect model produced errors.")
      warning("Linear mixed-effect regression produced an error. Verify the regression object of the result.\n")
      resultid <- NULL
      resultmean <- NULL
      regression <- model
      estimated <- NULL
    }

    #Renaming columns in $data object to original names
    intdata[, id := NULL]
    if(!is.null(resultid)){resultid[, id := NULL]}
    if(!is.null(estimated)){setnames(estimated, c("id_tmp", "time"), c(id, time))
      if(noinput){ #If there was no input, rename signal_estimated column in table estimated
        setnames(estimated,"signal_estimated",paste0(signal,"_estimated"))
      }}

    for(idx in seq(doreminames))
    {
      colnew <- doreminames[idx]
      colold <- originalnames[idx]
      tochange <- grep(colnew, names(intdata), value = T) #Composed names
      setnames(intdata, tochange, gsub(colnew, colold, tochange))
      if(!is.null(resultid)){
        if(colnew %in% c("id_tmp",doremiexc)){ #Only the id and the excitation columns need to change name
          tochange <- grep(colnew, names(resultid), value = T) #Composed names
          setnames(resultid, tochange, gsub(colnew,colold,tochange))
        }
      }
      if(!is.null(resultmean)){
        if(colnew %in% doremiexc){ #Only the excitation columns need to change name
          tochange <- grep(colnew, names(resultmean), value = T) #Composed names
          setnames(resultmean, tochange, gsub(colnew,colold,tochange))
        }
      }
    }
  }

  # Regression for SINGLE individuals ----------------------------------------
  else{ #If id is null it is assumed that all the data comes from a SINGLE INDIVIDUAL.
    print("Status: Time series data")
    #Saving original names and then renaming data table to internal names
    originalnames <- c(time, signal, input)
    doremiexc <- paste0("input", seq(input))
    doreminames <- c("time", "signal", doremiexc)
    setnames(intdata, originalnames, doreminames)

    #Calculation of the signal rollmean and first derivative of the signal column
    #Paste is only used to generate new column names based on the orignal ones (concatenate strings)
    intdata[, signal_rollmean := calculate.gold(signal, time, embedding)$dsignal[, 1]]
    intdata[, signal_derivate1 := calculate.gold(signal, time, embedding)$dsignal[, 2]]
    intdata[, time_derivate := calculate.gold(signal, time, embedding)$dtime]

    #Calculation of the roll mean of the excitation columns if there is at least one input column
    if (!noinput){
      myfun <- function(x){x[] <- c(rollmean(x, (embedding)), rep(NA, embedding - 1)); x}
      intdata[, (paste0(doremiexc,"_rollmean")) := lapply(.SD, myfun), .SDcols = doremiexc]
     }

    #In this case it is a linear regression and function lm is used instead of lmer
    #In the analysis for one individual there is no resultid table
    resultid <- NULL

    if(noinput){ # if there is no excitation signal
      model <- tryCatch({lm(signal_derivate1 ~ signal_rollmean, data = intdata)}, error = function(e) e)
      print("Status: Unknown excitation. Linear regression calculated")

    }
    else{ # if there is one or several excitation signals
      model <- tryCatch({lm(paste0("signal_derivate1 ~ signal_rollmean + ", paste(doremiexc, "rollmean ", collapse = "+", sep = "_")),
                             data = intdata)}, error = function(e) e)
      print("Status: One or several excitations. Linear regression calculated")
    }
    if (!inherits(model,"error")){ # if the regression worked
      print("Status: Linear mixed-effect model had no errors")
      summary <- summary(model) # Summary of the regression
      regression <- summary

      # Generate mean results
      # calculate the damping time: -1/damping_coeff
      resultmean <- setDT(list(dampingtime = -1L / summary$coefficients["signal_rollmean", "Estimate"]))

      # Extract the equilibrium value
      resultmean[, eqvalue := summary$coefficients["(Intercept)", "Estimate"] * resultmean[, dampingtime]]

      # Generation of the estimated signal using remi generate FOR ONE INDIVIDUAL with the estimated coefficients
      if (noinput){ #there is NO excitation as input
        print("Unknown excitation. Calculation of estimated signal through exponential fit")
        #Exponential model is assumed when no input is provided and thus a new fit is necessary:
        #log(y-B) = gamma*t + log(A) --> LINEAR EQUATION
        #B is known, it is the intercept*thao, from the model fit with the derivative inside the analysis function
        #A is unknown. It can be found by fitting log(y-B)~ t
        #We'll call the coefficients resulting from this new fit Ap and Bp. Ap=gamma, Bp=log(A)
        y <- intdata$signal
        t <- intdata$time
        B <- resultmean[, eqvalue] #intercept calculated previously with lm
        expmodel <- lm(log(abs(y-B)) ~ t)
        #Then it is assumed that the signal follows a decreasing exponential:
        # y = A* exp(gamma*t)+B
        if (!is.na(resultmean[, dampingtime]) && resultmean[, dampingtime] > 0){
          estimated <- setDT(list(time = intdata$time, signal_estimated =
                                  exp(expmodel$coefficients[1]) * exp(resultmean[, -1L/dampingtime] * time) + resultmean[, eqvalue]))
        }else{estimated <- NULL}

      }else{ #There is an excitation
        print("Status: One or several excitations. Calculation of estimated signal through convolution")
        intdata[, totalexc := 0]
        for (i in 1:length(input)) #For loop to go through all the inputs
        {
          resultmean[, c(paste0(doremiexc[i], "_coeff")) := summary$coefficients[paste0(doremiexc[i], "_rollmean"), "Estimate"] * resultmean[, dampingtime]]
          intdata[, totalexc:= totalexc + (summary$coefficients[paste0(doremiexc[i], "_rollmean"), "Estimate"]) * get(doremiexc[i])]
        }
        # Estimated table contains expanded time vector, minimum an maximum generated signals (for the two extreme scenarios of expanded excitation)
        # Expanded vectors according to deltat chosen
        # Generation of the expanded excitation vector according to desired deltat. Minimum and maximum values
        # Time vector and excitation vectors min and max
        deltat <- 0.1 * min(diff(intdata$time))
        #Expanded time vector: takes the minimum and the maximum of all the time intervals and creates the vector with deltat time intervals
        estimated <- intdata[, list(time = seq(floor(min(time, na.rm = T)), ceiling(max(time, na.rm = T)), deltat))]

        leftidx <- findInterval(intdata[, time], estimated[, time])

        #tmpsignal contains the value of the total excitation were the original times are found in the new time vector
        #and NA in the rest of the positions
        estimated[leftidx, tmpsignal := intdata[, totalexc]]


        #The na.locf function (last observation called forward) will repeat the last non NA value.
        #We use it to repeat the values in the excitation function to the left and to the right
        estimated[, exc_min := na.locf(tmpsignal, na.rm = F, fromLast = F)]
        estimated[, exc_max := na.locf(tmpsignal, na.rm = F, fromLast = T)]
        estimated <- estimated[!is.na(exc_min) & !is.na(exc_max)]

        #Removing tmpsignal from estimated and totalexc from intdata
        estimated[, tmpsignal := NULL]
        intdata[, totalexc := NULL]

        #Calculating the convolution
        estimated[, ymin := generate.remi(resultmean[, dampingtime], exc_min, time)$y + resultmean[, eqvalue]]
        estimated[, ymax := generate.remi(resultmean[, dampingtime], exc_max, time)$y + resultmean[, eqvalue]]
      }
    }else{ # if the regression didn't work, all coeffs will remain NA. The summary will contain the model and the estimated table wil be set to null
      print("Status: Linear regression produced errors.")
      warning("Linear regression produced an error. Verify the regression object of the result.\n")
      regression <- model
      resultmean <- NULL
      estimated <- NULL
    }
    #Renaming columns in $data object to original names
    if(!is.null(estimated)){#If there was no input, rename signal_estimated column in table estimated
      if(noinput){setnames(estimated,"signal_estimated",paste0(signal,"_estimated"))}
    }

    for(idx in seq(doreminames))
    {
      colnew <- doreminames[idx]
      colold <- originalnames[idx]
      tochange <- grep(colnew, names(intdata), value = T) #Composed names
      setnames(intdata, tochange, gsub(colnew, colold, tochange))
      if(!is.null(resultmean)){
        if(colnew %in% doremiexc){ #Only the excitation columns need to change name
          tochange <- grep(colnew, names(resultmean), value = T) #Composed names
          setnames(resultmean, tochange, gsub(colnew,colold,tochange))
        }
      }
    }
  }


  # Output the results of the function --------------------------------------
  #Excitation string
  #If there is no excitation term
  if (noinput){str_exc <- 0
  }else{str_exc <- input} # If there is one OR SEVERAL excitation columns

  res = list(data = intdata, resultid = resultid, resultmean = resultmean, regression = regression, estimated = estimated, embedding = embedding, str_time = time, str_exc = str_exc, str_signal = signal, str_id = id)
  class(res)= "doremi" # Class definition
  return(res)
}
# generate.panel.remi ----------------------------------------------

#' Simulation of various individual signals with intra and inter noise
#'
#' \code{generate.panel.remi} Generates signals with intra and inter individual noise for several individuals.

#' In order to do this, the function generates a pseudo-continuous signal per individual that is a solution to the first order differential equation:
#' \deqn{\frac{dy(t)}{dt} - \gamma y(t) = E(t)}
#' The analytical solution to this equation is a convolution between the Green function and the excitation term.
#' The function generates internally a pseudo-continuous signal to increase the precision with which the convolution is calculated. From this
#' expanded signal, the function samples points with a constant time step given by deltatf. These operations are repeated as many times as the value set in the input "nind". Once the signal is sampled,
#' intra-individual and inter-individual noise with normal distributions are added.

#' @param nind  number of individuals.
#' @param dampingtime Signal damping time. It corresponds to the time needed to reach
#' 37\% (1/e) of the difference between the equilibrium value and the amplitude of the signal reached when there is no excitation
#' (or 63\% of the maximum value for a constant excitation). It should be positive (dampingtime>0).
#' @inheritParams generate.excitation
#' @param internoise Is the inter-individual noise added. The dampingtime accross individuals follows a normal distribution centered on the input parameter dampingtime
#' with a standard deviation of internoise*dampingtime, except if any damping time is negative (see Details section).
#' @param intranoise Is the noise added in each individual signal. It also follows a normal distribution with a standard deviation equal to this parameter times the maximum
#' amplitude of the convolution when using an excitation containing a single pulse.

#' @keywords simulation, first order, differential equation
#' @return Returns a data frame with signal and time values starting at 0 and sampled at equal time steps deltatf in the time lapse tmax.
#' It contains the following columns:
#' \itemize{
#'    \item id - individual identifier (from 1 to nind).
#'    \item excitation - excitation signal generate through the generate.excitation
#'    \item time - time values
#'    \item dampedsignalraw - signal with no noise (inter noise added for each individual)
#'    \item dampedsignal - signal with intra noise added
#' }
#' @details Used for simulations in the context of the package.
#'
#' The function currently simulates only positive damping times corresponding to a regulated system. When the damping time is low
#' and the inter individual noise is high, some individuals' damping time could be negative. In that case, the damping time
#' distribution is truncated at 0.1*deltatf and values below are set to this limit. High values are symmetrically set at the upper percentile value
#' similar to a Winsorized mean. A warning provides the initial inter individual noise set as input argument and the inter individual
#' noise obtained after truncation.
#' @seealso \code{\link{generate.remi}} for calculation of the analytical solution to the differential equation.
#' Call the data frame $fulldata of the result for a full data frame with points generated at a very small deltatf in order
#' to build a pseudo-continuous function that will enhance the quality of the generated signal (see \code{\link{remi}}).
#' and \code{\link{generate.excitation}} for excitation signal generation
#' @examples
#' generate.panel.remi(nind = 5,
#'               dampingtime = 10,
#'               amplitude = c(5,10),
#'               nexc = 2,
#'               duration = 20,
#'               deltatf = 0.5,
#'               tmax = 200,
#'               minspacing = 0,
#'               internoise = 0.2,
#'               intranoise = 0.1)
#'@export
#'@importFrom data.table setDT
#'@importFrom data.table copy
#'@importFrom data.table :=
#'@importFrom data.table .N
#'@import stats
generate.panel.remi <- function(nind = 1,
                 dampingtime,
                 amplitude = 1,
                 nexc = 1,
                 duration = 10,
                 deltatf = 0.5,
                 tmax,
                 minspacing = 10,
                 internoise = 0,
                 intranoise = 0){
  #Internal time step to generate pseudo-continuous function
  deltat <- 0.1 * deltatf  #We do at least 10 divisions of the final time step chosen
  npoints <- tmax / deltat + 1

  # Generate simulation data for a given excitation and damping time
  # Generates id column (id being the individual number)
  # Creates a data table with id being the first column containing as many lines per individual as time points marked in npoints

  data <- setDT(list(id = unlist(lapply(c(1:nind), function(x){rep(x, npoints)}))))

  #Add to that data table an "excitation" column, containing the values of the excitation signal.
  #Creates a new excitation signal for each individual
  data[, excitation := generate.excitation(amplitude, nexc, duration, deltat, tmax, minspacing)$exc, by = id]
  data[, time := generate.excitation(amplitude, nexc, duration, deltat, tmax, minspacing)$t, by = id]

  #Creates a damping time vector by taking dampingtime input value and adding the internoise in a normal distribution
  #Contains as many elements as individuals
  dampingtimevec <- rnorm(nind, mean = dampingtime, sd = internoise * dampingtime)

  #If any value of the damping time vector is negative, the original value is used instead at that position of the vector
  if (any(dampingtimevec <= 0)){
    a <- 0.1 * deltatf #Threshold to truncate the dampingtime distribution
    perc <- as.vector(prop.table(table(dampingtimevec < a)))[2] #Calculation of the percentile of elements that are <0.1*deltatf (damping times of 0 are thus also excluded)
    b <- as.vector(quantile(dampingtimevec, probs = 1 - perc)) #Calculation of the symmetrical threshold

    #Truncating the normal distribution to these two thresholds
    dampingtimevec[dampingtimevec < a] <- a
    dampingtimevec[dampingtimevec > b] <- b

    #Calculating new sd and internoise of the truncated distribution
    nsd <- sd(dampingtimevec)
    ninternoise <- nsd / dampingtime

    warning("Some values for dampingtime where negative when adding internoise. Negative damping times imply signals
            that increase exponentially (diverge). This model generates signals that are self-regulated and thus,
            the normal distribution used has been truncated. The inter-individual noise added is of ", round(ninternoise,2), " instead of ",internoise,".\n")
  }

  #Creates the signals for each individual taking the damping time for that individual from dampingtimevec
  #dampedsignalraw is the signal WITHOUT NOISE
  data[, dampedsignalraw := generate.remi (dampingtimevec[.GRP], excitation, time)$y, by = id ]


  #Creates the signal for each individual with intra noise
  #In order to avoid adding an increased intra noise that will "grow" with each new pulse of the excitation signal, first
  #a "amplitudenorm" function is calculated in which the excitation has only one pulse that starts at the beginning of the time sequence
  #The intranoise is added regarding the maximum value of that function
  #As the excitation function can receive an amplitude vector and a duration vector, in order to normalize, it will
  #be necessary to pick the maximum of the amplitude and the excitation
  if (length(amplitude) > 1){amp <- max(amplitude)}
  else {amp <- amplitude}
  if (length(duration) > 1){dur <- max(duration) / deltat + 1} #As done in the excitation function
  else{dur <- duration / deltat + 1}

  data[, amplitudenorm := generate.remi(dampingtimevec[.GRP], rep(c(amp, 0), c(dur, (length(excitation)-dur))), time)$y, by = id ]
  data[, dampedsignal := dampedsignalraw + rnorm(.N, mean = 0, sd = intranoise * max(abs(amplitudenorm))), by = id ]

  #Returning fulldata
  fulldata <- copy(data)

  #Selecting equally spaced data with deltatf
  #definition of time step for equally spaced values.
  data <- data[time %in% seq(0, tmax, deltatf), .SD,by = id] #Keeping only values from data table that are in the time intervals chosen

  #Defining components of the result value
  res <- list(fulldata = fulldata[, !"amplitudenorm", with = FALSE], data = data[, !c("amplitudenorm"), with = FALSE])
  class(res)= "doremidata" #Class definition
  return(res)
}
