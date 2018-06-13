# calculate_gold ----------------------------------------------------------
#' Calculation of derivatives through using the Gold method
#'
#' \code{Calculate_gold} calculates the derivatives of a function according to the Gold method described in \url{https://doi.org/10.1080/00273171.2010.498294}
#' including management of variable time steps.
#' @param signal is a vector containing the values of the function to be derivated.
#' @param time Is a vector containing time values.  It has the same length as signal.
#' @param embedding is a signed integer indicating the number of points to consider for derivative calculation. Embedding must be greater than 2 because at least
#' two points are needed for the calculation of the first derivative and at least 3 for the calulation of the second derivative. If it is not the case, the
#' function displays an error message. If not specified, embedding=2 by default.
#' @keywords derivative, embed, rollmean
#' @return Returns a list containing three columns:
#' dtime- contains the time values in which the derivative was calculated. That is, the roll mean of the input "time" vector in "embedding" points.
#' dsignal- contains the calculated derivative of the "signal" vector according to the GOLD method, considering "embedding" points.
#' embedding- contains the number of points used for the derivative calculation and time roll means, which is constant.
#'
#' @examples
#' time <- c(1:500)/100
#' signal <- time^2
#' result <- calculate_gold(signal = signal, time = time, embedding = 5)
#'@export
#'@importFrom zoo rollmean
calculate_gold <-  function(signal,
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
    #as groups of time points the embed funcion gives (thus, the number of lines
    #of the tembed matrix) Columns: It will have 3 columns because: First column
    #will contain the signal mean value in the time points considered in
    #"embedding" Second column will contain the signal derivative in those same
    #time points Third column will contain the signal second derivative in those
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
    # Addition of NA so that the derivative rows are the same as those of timeserie

  } else if (embedding == 2){
    warning("Only first derivative can be calculated with an embedding of 2.\n")
    derivative <- cbind(rollmean(signal, embedding), diff(signal) / diff(time))
    #Appends the two columns: 1. mean time values and 2. The span calculated as
    #the difference of two signal values (going forward) divided by the time
    #interval
    derivative <- rbind(derivative,matrix(data = NA, ncol = 2, nrow = embedding - 1))
    # Addition of NA so that the derivative rows are the same as those of signal
  } else{
    stop("Embedding>=2 for the calculation of a first derivative and 3 for a second derivative.\n")
  }
  time_derivative <-c(rollmean(time, embedding), rep(NA, embedding - 1))
  # Addition of NA so that the time rows are the same as those of signal

  returnobject <- list("dtime" = time_derivative,
                       "dsignal" = derivative,
                       "embedding" = embedding)
  return(returnobject)
}


# excitation_function -----------------------------------------------------
#' Excitation signal generation
#'
#' \code{excitation_function} returns a vector containing the values of a generated pulsed excitation signal with the pulses randomnly placed
#' for a given amplitude, duration and spacing between the pulses.
#' @param amplitude  is a signed integer or vector of signed integers indicating the amplitude of the pulse. It must be greater than 0,
#' otherwise an error is generated. It should also have as many elements as the number of pulses (nexc). If this is not the case and the
#' elements are less than the number of pulses, the amplitude vector will be "recycled" and the elements from it will be repeated until
#' all the pulses are covered (for instance, if the signal has 6 pulses and the amplitude vector has two elements, pulses 1,3 and 5 will
#' have the same amplitude as the first element of the amplitude vector and pulses 2,4 and 6 that of the second element).A warning will also
#' be generated.
#' @param nexc is a signed integer indicating the number of pulses to be generated. It must be greater than 0, otherwise an error is generated.
#' @param duration is a signed integer or vector of signed integers indicating the duration of the pulse in time units. It must be greater than 0,
#' otherwise an error is generated. As for the amplitude, it should also have as many elements as the number of pulses (nexc). If this is not the
#' case and the elements are less than the number of pulses, the amplitude vector will be "recycled" and the elements from it will be repeated until
#' all the pulses are covered. A warning will also be generated.
#' @param deltatf is a signed decimal indicating the time step.
#' @param tmax is a signed integer indicating the maximum time range of the excitation in time units. The time vector generated will go from 0 to tmax.
#' @param minspacing as pulses are generated randomnly, minspacing is a signed integer that indicates minimum spacing between pulses, thus avoiding
#' overlapping of the pulses in time. It can be 0 indicating that pulses can follow one another.
#' @keywords excitation, simulation
#' @return Returns two vectors:
#' E- vector containing the values of the excitation generated.
#' t- vector containing the values of time generated.
#' @details Used for simulations in the context of the package. Beware that the following condition should apply:
#' \deqn{tmax >= (duration+minspacing)*nexc}
#' so that the pulses "fit" in the time lapse defined.
#' The package "seewave" already contains a function that generates a pulse (see function pulsew). The improvement of this function, regarding pulsew is that
#' it is possible to generate pulses of different duration and amplitude, and not just unitary (amplitude=1).
#'
#' @examples
#' excitation_function(c(1,10,20),3,c(1,2,4),0.5,100,10)
#' excitation_function(3,6,2,1,200,2)
#'@export
#'@importFrom utils head
# Excitation signal generation.
excitation_function <- function(amplitude = 1,
                                nexc = 1,
                                duration = 10,
                                deltatf = 1,
                                tmax = 100,
                                minspacing = 10){
  #Error management
  if (any(duration <= 0) | any(amplitude == 0) | nexc <= 0){
    stop("Invalid input parameters. At least one excitation must be defined. Duration and nexc must be greater than 0.Amplitude must be different from 0.\n")
  }
  if (nexc < length(duration) | nexc < length(amplitude)){
    warning("The number of excitations nexc is smaller than the number of elements in amplitude and duration. Only the first elements of these vectors were considered.\n")
  }

  if (nexc > length(duration) && length(duration) > 1){
    duration <- rep(duration, ceiling(nexc / length(duration)))
    duration <- duration[1:nexc]
    warning("The number of excitations nexc was higher than the durations defined. The values from the vector were repeated.\n")
  }
  if (nexc > length(amplitude) && length(amplitude) > 1){
    amplitude <- rep(amplitude, ceiling(nexc / length(amplitude)));
    amplitude <- amplitude[1:nexc]
    warning("The number of excitations nexc was higher than the amplitudes defined. The values from the vector were repeated.\n")
  }
  if (tmax < (sum(duration) + minspacing) * nexc){
    stop("Non valid parameters. tmax should be greater than (duration+minspacing)*nexc.\n")
  }

  if (length(duration) == 1) {dur <- duration * nexc}
  #If it's a scalar, it assigns the value directly. If it is not, it takes the value from the vector
  #at the position of the correponding excitation
  else {dur <- sum(duration)}
  trest <- tmax - dur - minspacing * (nexc - 1) #at the beginning, calculates the time left (to calculate random time from it) is tmax-total pulse duration - total minspacing duration
  cumt <- 0 #cumulated time initialized to 0

  for (i in 1:nexc){
    #Generation of a single unitary pulse (single pulse with amplitude=1)
    if (length(duration) == 1) {dur <- duration; durleft <- dur * (nexc - i)} #If it's a scalar, it assigns the value directly. If it is not, it takes the value from the vector
    #at the position of the correponding excitation
    else {dur <- duration[i]; durleft <- sum(duration[(i+1):nexc])}

    tal<-tmax #Random time is initialized so that it enters the while loop the first time

    while (tal > trest){ #While the random time coming from the sample is greater than the time left before tmax
      if (i == 1){tal <- sample(0:trest, 1, replace = T)}
      else{
        if(trest == minspacing){tal <- minspacing}
        else{tal <- sample(minspacing:trest, 1, replace = T)}
      }
    }
    nf<-tmax/deltatf+1
    sp <- rep(c(0, 1), c(tal / deltatf + 1, dur / deltatf + 1))
    #simple unitary pulse
    cumt <- cumt + tal + dur #cumulated time

    #Generation of amplified simple pulse
    if (length(amplitude) == 1) {amp <- amplitude} #Same as for duration.
    else {amp <- amplitude[i]}
    sp <- sp * amp

    #Append calculated single pulse to result vector
    if (i == 1){E <- sp}
    else{E <- append(E, sp)}

    #Calculating the time left so that a new sample of random time can be taken
    trest <- tmax - cumt - durleft - minspacing * (nexc - i - 1)

  }

  #Verify E vector length
  if (length(E) < nf) {
    E <- append(E, rep(0, nf -length(E))) #If the length is smaller than npoints, fill with 0
  }
  else{
    E <- E[1:nf]
    warning("Due to input parameters introduced, vector size was larger than npoints and was cut to this value.\n")
  }
  #Generation of time vector
  tim <- seq(0, tmax, deltatf)

  data <- list(y = E, t = tim)
  return(data)
}



# doremi_analyse_order1 -----------------------------------------------------


#' DOREMI first order analysis function
#'
#' \code{doremi_analyse_order1}  the function aims to quantify how close is a study signal to the solution of the first order linear
#' constant coefficient differential equation of the form: \deqn{\dot{y}(t) = -\gamma y(t) + E(t)}
#' Where y(t) is the individual's signal and E(t) is the excitation.
#' In order to do this, the function uses multilevel regression and calculates the coefficients, fixed and random effecs and the estimated signal.
#' @param userdata Is a data frame or data table containing at least two columns: one containing the identifier for each individual and the other one the
#' signal to be analysed.
#' @param id Is a STRING containing the NAME of the column of userdata containing the idENTIFIER of the individual.
#' If this parameter is not entered when calling the function, a single individual is assumed and a linear regression is done instead.
#' @param input Is a STRING or a VECTOR OF STRINGS containing the NAME of userdata column/s containing the EXCITATION vector/s.
#' If this parameter is not entered when calling the function,
#'         the excitation is assumed to be constant/non existent. In this case, the multi-level regression is still carried out but no coefficient is calculated
#'         for the excitation term. The function then uses the parameters estimated by the regression to carry out an exponential fit of the signal.
#' If this parameter is a vector, the function will consider an excitation column each column of userdata whose name is contained in the input vector.
#'         The function will then include all the excitations in the multilevel regression and calculate a coefficient for each one of these.
#' @param time  Is a STRING containing the NAME of the column of userdata containing the time vector. If this parameter is not entered when calling the function,
#' it is assumed that time steps are of 1 unit and the time vector is generated internally in the function.
#' @param signalcolumn Is a STRING containing the NAME of the column of the data frame containing the SIGNAL to be studied.This parameter can't be left empty.
#' @param embedding Is a positive integer containing the number of points to be used for the calculation of the derivatives. Its value by default is 2 as at
#' least two points are needed for the calculation of the first derivative.
#' @keywords analysis, first order, exponential
#' @return The function returns five objects:
#' \enumerate{
#'  \item data- A data.frame including the input data (for comparison purposes) and other intermediate calculations used to prepare the variables for
#'   the fit, more specifically:
#'     dampedsignal_rollmean - calculation of the moving average of the signal in embedding points.
#'     dampedsignal_derivate1 - calculation of the first derivative of the signal according to the GOLD method in embedding points.
#'     timecol_derivate - calculation of the moving average of the time vector in embedding points.
#'     excitation_rollmean - calculation of the moving average of the excitation vector in embedding points.
#'  \item resultid- A data.frame including for each individual, listed by his/her id number, the damping time, the excitation coefficient and the
#'  equilibrium value (see variables presented above).
#'  \item resultmean- A data.frame including the average of the fixed effects of the three coefficients mentioned above for all the individiuals.
#'  \item regression- A list containing the summary of the multilevel regression (random and fixed parts, residuals, t values).
#'  \item estimated- A data.frame containing the estimated signal calculated as the convolution of the Green function and the excitation
#'  vector with an added offset (see above). There are two extreme cases in the generation of the signal and these depend on sampling. The excitation
#'  vector is expanded to generate a pseudo-continuous signal and increase accuracy when calculating the convolution. Missing data in the exitation signal
#'  is thus completed taking the previous known value (to the "left" and thus generating the signal exc_min in the data.frame) or taking the next known value
#'  (to the "left" and thus generating the signal exc_max in the data.frame). With these two expanded excitation vectors, the two exreme cases of the
#'  estimated signal are calculated by carrying out the convolution between the Green function (decreasing exponential with the damping time calculated by
#'  the multilevel regression) and each one of the excitation signals described, thus generating ymin and ymax respectively.
#'}

#' @details The analysis performs the following multilevel regression:
#' \deqn{y_{ij}'  \sim   b_{0} +b_{0i}+b_{1} y_{ij}+b_{2} E_{ij}+u_{1j} y_{ij}+u_{2j} E_{ij}+e_{ij}}
#' with:
#' j accounting for the different individuals
#' i for the time
#' \eqn{e_{ij}} for the error term
#' \eqn{y_{ij}'} is the derivative calculated on \eqn{N_{E}} points
#' y and E are the signal and the excitation averaged on \eqn{N_{E}} points.
#' From this equation, the coefficients provided as outputs of the function and that characterize the signal are calculated as follows:
#' \itemize{
#'   \item Damping time:  \eqn{\tau _{i} =  \frac{1}{ \gamma _{i} }}  with \eqn{\gamma _{i} =  b_{1} + u_{1j} }
#'   \item Excitation coefficient: \eqn{\frac{b_{2} + u_{2j}}{\gamma _{i}}}It is the proportionality between the excitation and the maximum value reached
#' by the variable experiencing this excitation.
#'   \item Equilibrium value: \eqn{\frac{b_{0} + b_{0i} }{\gamma _{i}}} It is the stable value reached in the absence of excitation.
#' }
#' The estimation is performed using the funtion lmer if there are several indiviuals or lm if there is only one.
#' With the above estimated parameters, the estimated signal can be reconstructed for
#' each individual by first performing the convolution of the excitation with the Green function , having a damping rate equal to the estimated
#' and then offsetting the resulting signal with the equilibrium value.
#'
#' DOREMI can be used to analyse several situations in different domains, for instance:
#' \itemize{
#'   \item \strong{Medicine}
#'          Excitation: load in a resistive bicycle during an effort test
#'          Signal: individual's cardiac rythm
#'   \item \strong{Sociology}
#'          Excitation: occurrence of perturbing events in an individual's life
#'          Signal: measurement of biomarkers in the individual with respect to time
#'   \item \strong{Psychology}
#'          Excitation: Sessions and trials where the individuals are required to carry out a task
#'          Signal: response time for the execution of the task
#' }
#' @seealso \code{\link{simulation_generate_order1}} to generate simulation signals that satisfy the equation presented
#' @examples
#' mydata <- simulation_generate_order1(nindividuals = 5,
#'                            dampingtime = 10,
#'                            amplitude = c(5,10),
#'                            nexc = 2,
#'                            duration = 20,
#'                            deltatf = 2,
#'                            tmax = 200,
#'                            minspacing = 0,
#'                            internoise = 0.2,
#'                            intranoise = 0.1)
#' myresult <- doremi_analyse_order1(userdata = mydata$data,
#'                            id = "id",
#'                            input = "excitation",
#'                            time = "timecol",
#'                            signalcolumn = "dampedsignal",
#'                            embedding = 5)
#'@export
#'@importFrom data.table setDT
#'@importFrom data.table copy
#'@importFrom data.table setkey
#'@importFrom data.table shift
#'@importFrom data.table :=
#'@importFrom data.table .SD
#'@importFrom data.table .N
#'@importFrom data.table .GRP
#'@importFrom data.table setnames
#'@importFrom data.table setcolorder
#'@importFrom lme4 lmer
#'@importFrom lme4 lmerControl
#'@importFrom lme4 ranef
#'@importFrom stats embed
#'@importFrom stats lm
#'@importFrom zoo rollmean
doremi_analyse_order1 <- function(userdata,
                                  id = NULL,
                                  input = NULL,
                                  time = NULL,
                                  signalcolumn,
                                  embedding = 2){
  data <- copy(userdata) #Keeps the original of data so that it can rename columns freely
  data <- setDT(data) # Convert input data to data.table.
  noinput <- FALSE #Flag that will allow to differentiate if there is an excitation term or not when doing the regression
  #Error management
  #Verifying column names repeated in data table.
  if(any(duplicated(colnames(data)))){stop("input datatable contains duplicated column names. Please correct these in order to launch the analysis function.\n")}

  #If id,input,time,signalcolumn are not strings containing the name of columns in "data", stop the function
  if(!is.null(id) && !is.character(id)){stop("id should be a string containing the name of the column in data that contains the individual identifier.\n")}
  if(!is.null(input) && !is.character(input)){stop("input should be a string containing the name of the column in data that contains the excitation.\n")}
  if(!is.null(time) && !is.character(time)){stop("time should be a string containing the name of the column in data that contains the time.\n")}
  if(!is.character(signalcolumn)){stop("signalcolumn should be a string containing the name of the column in data that contains the signal of the individual.\n")}

  #Verifying if the excitation is given as input. If not, it is set to 0. In this case, a warning will be generated
  #Indicating that the input has been set to 0.
  # It is difficult to verify if the signal is constant all along as it can be constant during the pulse
  if (is.null(input)){data[, inputcol := 0]
    input = "inputcol"
    noinput<-TRUE #This flag will be needed later as if there is no input, coefficients for the excitation term will not be calculated in the regression
    warning("No excitation signal introduced as input. input was set to 0.\n")
  }

  if (is.null(time)){data[, timecol := 1L * c(1:.N), by = id]
    time <- "timecol" # if no time set it to a 1 sec step vector
    warning("No time vector introduced as input. A 1 unit increment time vector is generated.\n")
  }

  #Converting data if columns are of type "factor" or "string" instead of numeric
  if (any(is.factor(sapply(data, class))) | any(is.character(sapply(data, class)))){
    strtmp <- c(input, time, signalcolumn)
    data[, strtmp := lapply(.SD, function(x) {as.numeric(as.character(x))}), .SDcols = strtmp]
    warning("Some columns were found to be of the factor/string type and were converted to numeric.\n")
  }

  #Rename column containing id to "id" because further functions use this identifier to work. Also display warning message saying that column has been renamed.
  if (!is.null(id) && id!="id"){

    #Setting some of the column names (necessary for data treatment)
    setnames(data, old = id, new = "id")

  }

  if (embedding == 2){warning("Only first derivative can be calculated with an embedding of 2.\n")}

  #Suppress NA elements if there is an NA in time or in signalcolumn
  #This suppresses the entire row of the data table
  data <- data[!is.na(get(time)) & !is.na(get(signalcolumn))]

  #If in the data left, there are some NA values in the excitation vector, set them to zero.
  #This is to avoid loosing data from signal and time vectors, as sometimes when people fill in the signal data they put NA to mean no excitation, thus 0.
  na.to.0 <- function(x){x[is.na(x)] <- 0; x}
  data[, (input) := lapply(.SD, na.to.0), .SDcols = input]

  #Calculation of the signal rollmean and first derivative of the signal column
  #Paste is only used to generate new column names based on the orignal ones (concatenate strings)
  options(warn = -1) #Turns off warnings to avoid the calculate_gold function activating warnings every time it is called

  data[, c(paste0(signalcolumn,"_rollmean")) := calculate_gold(get(signalcolumn),get(time),embedding)$dsignal[,1], by = id]
  data[, c(paste0(signalcolumn,"_derivate1")) := calculate_gold(get(signalcolumn),get(time),embedding)$dsignal[,2], by = id]
  data[, c(paste0(time,"_derivate")) := calculate_gold(get(signalcolumn),get(time),embedding)$dtime, by = id]


  options(warn = 0) #Turns warnings back on

  #Calculation of the roll mean of the excitation columns if there is at least one input column
  if (!noinput){
    myfun <- function(x){x[] <- c(rollmean(x, (embedding)), rep(NA,embedding - 1)); x}
    data[, (paste0(input,"_rolled")) := lapply(.SD, myfun), .SDcols = input, by = id]

    #Calculation of the final excitation column (superposition of all the excitation signals)
    data[, excitation_sum := unlist(Reduce(function(a,b) Map(`+`,a,b), .SD)), .SDcols = input]
  }

  #Result data frames
  #The first table contains the input data
  #The second table contains the mean results for gamma and thao for each individual
  resultid <- setDT(list( id = unique(data$id) ))
  #The third table contains the mean values for gamma and thao for all the individuals (single line)
  #Generate mean results with convergence criterions
  resultmean <- setDT(list(id = "All"))

  # Regression for SEVERAL INDIVidUALS --------------------------------------


  # First order derivative equation mixed regression
  if (!is.null(id)){ #If the id column is not null, it is assumed that the data comes from SEVERAL INDIVidUALS
    #whose identifiers are in the id column.

    setkey(resultid, id) #sorts the data table by id

    if (noinput){ # if there is no excitation signal
      model <- tryCatch({lmer(paste0(signalcolumn,"_derivate1 ~ ", signalcolumn, "_rollmean + (1 + ", signalcolumn, "_rollmean |id)"),
                              data = data, REML=TRUE,
                              control = lmerControl(calc.derivs = FALSE,optimizer = "nloptwrap"))}, error = function(e) e)
    }else{ # if there is one OR SEVERAL excitation signals
      model <- tryCatch({lmer(paste0(signalcolumn, "_derivate1 ~ ", signalcolumn, "_rollmean + (1 + ", paste(input, "rolled ", collapse = "+", sep = "_")," + ",
                                     signalcolumn, "_rollmean |id) + ", paste(input, "rolled ", collapse = "+",sep = "_")),
                              data = data, REML=TRUE,
                              control = lmerControl(calc.derivs = FALSE, optimizer = "nloptwrap"))}, error = function(e) e)

    }
    if (!inherits(model,"error")){ # if the regression worked
      summary <- summary(model) # Summary of the regression
      random <- ranef(model) # Variation of the estimate coefficient over the individuals
      regression <- list(summary, random) # list to output both results: summary, and the table from ranef

      # calculate the damping time for all signal columns
      resultmean[, c(paste0(signalcolumn, "_dampingtime")) :=
                   -1L/summary$coefficients[paste0(signalcolumn, "_rollmean"), "Estimate"]] # calculate the damping time: -1/damping_coeff

      # Extract the intercept coeff (equilibrium value)
      resultmean[, c(paste0(signalcolumn,"_eqvalue")) := summary$coefficients["(Intercept)","Estimate"] * resultmean[, get(paste0(signalcolumn, "_dampingtime"))]]


      #Generate the results for each individual (second output table)
      # Calculate the damping time from the damping coeff

      # Damping time in resultid ------------------------------------------------


      resultid[, c(paste0(signalcolumn,"_dampingtime")) := -1L/(summary$coefficients[paste0(signalcolumn,"_rollmean"),"Estimate"] + random$id[.GRP,paste0(signalcolumn,"_rollmean")]), by = id]

      # Extract the intercept (equilibrium value) calculated for each individual (present in random, regression table)
      # Offset in resultid ------------------------------------------------------
      resultid[, c(paste0(signalcolumn,"_eqvalue")) := (summary$coefficients["(Intercept)", "Estimate"] + random$id[.GRP, "(Intercept)"])*resultid[.GRP, get(paste0(signalcolumn, "_dampingtime"))], by = id]

      # Generation of the fitted signal for all id using remi generate FOR SEVERAL INDIVidUALS


      if(noinput){
        #Exponential model is assumed when no input is provided and thus a new fit is necessary BY INDIVidUAL
        #log(y-B) = gamma*t + log(A) --> LINEAR EQUATION
        #B is known, it is the intercept*thao, from the model fit with the derivative inside the analysis function
        #A is unknown. It can be found by fitting log(y-B)~ t
        #We'll call the coefficients reulting from this new fit Ap and Bp
        #Ap=gamma
        #Bp=log(A)

        data[, expfit_A := lm(log(get(signalcolumn) - resultid[.GRP, get(paste0(signalcolumn, "_eqvalue"))]) ~ get(time))$coefficients[1], by = id]
      }
      else{ # If there is an excitation term
        # Extract the excitation coeff for each excitation
        data[, totalexc := 0]
        for (i in 1:length(input)){ #For loop to go through all the inputs
          resultmean[, c(paste0(input[i],"_exccoeff")) := summary$coefficients[paste0(input[i], "_rolled"), "Estimate"] * resultmean[, get(paste0(signalcolumn, "_dampingtime"))]]


          #If variation of the excitation coefficient accross individuals needed:
          #And for each individual: the mean coeff (sumary$coeff) + the variation per Individual (in random)

          # Excitation coefficient in resultid --------------------------------------

          resultid[, c(paste0(input[i],"_exccoeff")) :=
                     (summary$coefficients[paste0(input[i], "_rolled"), "Estimate"] + random$id[.GRP,paste0(input[i], "_rolled")]) * resultid[.GRP, get(paste0(signalcolumn, "_dampingtime"))], by = id]
          data[, totalexc := totalexc+(summary$coefficients[paste0(input[i], "_rolled"), "Estimate"] + random$id[.GRP,paste0(input[i], "_rolled")]) * get(input[i]), by = id]
        }
      }

      # Write any error message from the regression
      resultmean[, c(paste0(signalcolumn,"_fitmsg")) := summary$fitMsgs]

      # Write a warning if any of the damping times calculated was negative
      if (any(is.na(resultid[, get(paste0(signalcolumn,"_dampingtime"))])) | any(resultid[, get(paste0(signalcolumn,"_dampingtime"))] < 0)){
        warning("Some of the damping times calculated were negative and thus, the estimated signal is not generated for these.\n")
      }
      if (noinput){ #There is no excitation signal as input
        data[, c(paste0(signalcolumn,"_estimated")) :=
               if(!is.na(resultid[.GRP, get(paste0(signalcolumn,"_dampingtime"))]) && resultid[.GRP, get(paste0(signalcolumn,"_dampingtime"))] > 0){
                 # if there is a damping time that has been calculated and if it is greater than 0 (decreasing exponential)
                 #and if there is an excitation coefficient (excitation term was found in the inputs)

                 #Then it is assumed that the signal follows a decreasing exponential:
                 # y = A* exp(gamma*t)+B
                 #expmodel comes from fitting log(y-B)~t. Calculated above
                 exp(expfit_A) * exp(-1L / resultid[.GRP, get(paste0(signalcolumn,"_dampingtime"))] * get(time)) + resultid[.GRP, get(paste0(signalcolumn,"_eqvalue"))]

               }else{NaN}, by = id]
        estimated <- NULL
      }else{
        # Generation of the third result table called \"estimated\"
        # Contains expanded time vector, minimum an maximum generated signals (for the two extreme scenarios of expanded excitation)

        #Expanded vectors according to precision
        #Generation of the expanded excitation vector according to desired deltat. Minimum and maximum values

        precision <- 100 #number of divisions in each time step
        estimated <- data[, list(exc_min = rep(totalexc[1:(.N-1)], each = precision[1]), tmax = max(get(time)), intervals = diff(get(time))/precision), by = id]
        estimated[, timecol := c(0 + cumsum(c(0, intervals[1:(length(intervals)-1)]))), by = id]
        #Shifting exc_min precision-1 values to the LEFT to obtain exc_max (the parameter "lead" means shift to the left whereas "lag" means shift to the right)
        estimated[, exc_max := shift(exc_min, n = precision - 1, fill = 0, type = "lead"), by = id]

        #Generation of the estimated signals
        estimated[, ymin := doremi_generate_order1(resultid[.GRP, get(paste0(signalcolumn,"_dampingtime"))],exc_min,timecol)$y+
                    resultid[.GRP, get(paste0(signalcolumn,"_eqvalue"))], by = id]

        estimated[, ymax := doremi_generate_order1(resultid[.GRP, get(paste0(signalcolumn,"_dampingtime"))],exc_max,timecol)$y+
                    resultid[.GRP, get(paste0(signalcolumn,"_eqvalue"))], by = id]

        #Time vector and excitation vectors min and max
        deltat <- 0.1
        #Expanded time vector: takes the minimum and the maximum of all the time intervals and creates the vector with deltat time intervals
        estimated2 <- data[, list(timecol = seq(floor(min(get(time), na.rm = T)), ceiling(max(get(time), na.rm = T)), deltat)), by = id]

        #Finding the values of the original time vector in the expanded time vector by using the function "findInterval"
        IDvec <- unique(data$id)

        for (idx in seq(IDvec)){ #It is necessary to do a loop because findInterval finds the index in which the value is found in the original time vector
          #And it will be necessary to shift these indexes according to in which id we are
          leftidx <- findInterval(data[id == IDvec[idx], get(time)], estimated2[id == IDvec[idx], timecol])
          #tmpsignal contains the vlue of exctotal were the original times are found in the new time vector
          #and NA in the rest of the positions
          estimated2[nrow(estimated2[id < IDvec[idx]]) + leftidx, tmpsignal := data[id == IDvec[idx], totalexc]]
       }

        #The na.locf function (last observation called forward) will repeat the last non NA value.
        #We use it to repeat the values in the excitation function to the left and to the right
        estimated2[, exc_min := na.locf(tmpsignal, na.rm = F, fromLast = F), by = id]
        estimated2[, exc_max := na.locf(tmpsignal, na.rm = F, fromLast = T), by = id]

        #Calculating the convolution
        estimated2[, ymin := doremi_generate_order1(resultid[.GRP, get(paste0(signalcolumn,"_dampingtime"))],exc_min,timecol)$y+
                     resultid[.GRP, get(paste0(signalcolumn,"_eqvalue"))], by = id]
        estimated2[, ymax := doremi_generate_order1(resultid[.GRP, get(paste0(signalcolumn,"_dampingtime"))],exc_max,timecol)$y+
                     resultid[.GRP, get(paste0(signalcolumn,"_eqvalue"))], by = id]

        estimated <- estimated[, !c("tmax","intervals")]
        setcolorder(estimated, c("id", "timecol", "exc_min", "exc_max", "ymin", "ymax"))
      }
    }else{ # if the regression didn't work, set to NA all coeffs
      resultid[, c(paste0(signalcolumn,c("_dampingtime","_exccoeff","_eqvalue"))) := NA]
      resultmean[, c(paste0(signalcolumn,c("_dampingtime","_exccoeff","_eqvalue"))) := NA]
      regression <- model
      estimated <- NULL
    }

    # Output the results for the function
    return(list(data = data[,!"strtmp"], resultid = resultid, resultmean = resultmean, regression = regression, estimated = estimated, estimated2 = estimated2))

  }

  # Regression for SINGLE INDIVidUAL ----------------------------------------


  else{ #If id is null it is assumed that all the data comes from a SINGLE INDIVidUAL.
    #In this case it is a linear regression and function lm is used instead of lmer
    if(noinput){ # if there is no excitation signal

      model <- tryCatch({lm(paste0( signalcolumn,"_derivate1 ~ ",signalcolumn,"_rollmean"),
                            data=data)}, error=function(e) e)

    }
    else{ # if there is one or several excitation signals
      model <- tryCatch({ lm(paste0( signalcolumn,"_derivate1 ~ ",signalcolumn,"_rollmean + ",paste(input,"rolled ",collapse="+",sep="_")),
                             data=data)}, error=function(e) e)

    }
    if (!inherits(model,"error")){ # if the regression worked
      summary <- summary(model) # Summary of the regression

      #Generate mean results (third output table, there will not be a second output table as there is a single individual)
      #with convergence criterions
      # Extract the damping coefficient from the regression summary
      resultmean <- setDT(list(id = "All"))

      # calculate the damping time
      resultmean[, c(paste0(signalcolumn, "_dampingtime")) :=
                   -1L / summary$coefficients[paste0(signalcolumn, "_rollmean"), "Estimate"] ] # calculate the damping time: -1/damping_coeff

      # Extract the intercept coeff (equilibrium value)
      resultmean[, c(paste0(signalcolumn, "_eqvalue")) := summary$coefficients["(Intercept)", "Estimate"] * resultmean[, get(paste0(signalcolumn, "_dampingtime"))]]

      # Extract the excitation coeff for each excitation
      if (!noinput){ # If there is an excitation term
        data[, totalexc := 0]
        for (i in 1:length(input)) #For loop to go through all the inputs
        {
          resultmean[, c(paste0(input[i], "_exccoeff")) := summary$coefficients[paste0(input[i], "_rolled"), "Estimate"] * resultmean[, get(paste0(signalcolumn, "_dampingtime"))]]
          data[, totalexc:= totalexc + (summary$coefficients[paste0(input[i], "_rolled"), "Estimate"]) * get(input[i]), by = id]
        }
      }


      # Write any error message from the regression
      resultmean[, c(paste0(signalcolumn, "_fitmsg")) := summary$fitMsgs ]

      # Generation of the fitted signal using remi generate FOR ONE INDIVidUAL with the estimated damping time

      if (noinput){
        #Exponential model is assumed when no input is provided and thus a new fit is necessary:
        #log(y-B) = gamma*t + log(A) --> LINEAR EQUATION
        #B is known, it is the intercept*thao, from the model fit with the derivative inside the analysis function
        #A is unknown. It can be found by fitting log(y-B)~ t
        #We'll call the coefficients reulting from this new fit Ap and Bp
        #Ap=gamma
        #Bp=log(A)

        y <- data[[signalcolumn]] #signal
        t <- data[[time]]
        B <- resultmean[, get(paste0(signalcolumn,"_eqvalue"))] #intercept calculated previously wih lm

        expmodel <- lm(log(y-B) ~ t)
      }
      if (noinput){ #there is NO excitation as input
        data[, c(paste0(signalcolumn, "_estimated")) :=
               if (!is.na(resultmean[, get(paste0(signalcolumn,"_dampingtime"))]) && resultmean[, get(paste0(signalcolumn,"_dampingtime"))] > 0 ){
                 #Then it is assumed that the signal follows a decreasing exponential:
                 # y = A* exp(gamma*t)+B
                 #expmodel comes from fitting log(y-B)~t. Calculated above
                 exp(expmodel$coefficients[1]) * exp(resultmean[, -1L/get(paste0(signalcolumn, "_dampingtime"))] * get(time)) +
                   resultmean[, get(paste0(signalcolumn, "_eqvalue"))]
               }else {NaN}, by = id]
        estimated <- NULL
      }else{ #There is an excitation

        # Generation of the third result table called \"estimated\"
        # Contains expanded time vector, minimum an maximum generated signals (for the two extreme scenarios of expanded excitation)

        #Expanded vectors according to precision
        #Generation of the expanded excitation vector according to precision. Minimum and maximum values
        precision <- 100 #number of divisions in each time step
        estimated <- data[, list(exc_min = rep(totalexc[1:(.N-1)], each = precision[1]), tmax = max(timecol), intervals = diff(timecol) / precision), by = id]
        estimated[, timecol := c(0 + cumsum(c(0, intervals[1:(length(intervals)-1)]))), by = id]

        #Shifting exc_min precision-1 values to the LEFT to obtain exc_max (the parameter "lead" means shift to the left whereas "lag" means shift to the right)
        estimated[, exc_max := shift(exc_min, n = precision - 1, fill = 0, type ="lead"), by = id]

        #Generation of the estimated signals
        estimated[, ymin := doremi_generate_order1(resultmean[, get(paste0(signalcolumn, "_dampingtime"))], exc_min, timecol)$y +
                    resultmean[, get(paste0(signalcolumn,"_eqvalue"))], by = id]

        estimated[, ymax := doremi_generate_order1(resultmean[, get(paste0(signalcolumn, "_dampingtime"))], exc_max, timecol)$y +
                    resultmean[, get(paste0(signalcolumn,"_eqvalue"))], by = id]
        estimated <- estimated[, !c("tmax", "intervals")]
        setcolorder(estimated, c("timecol", "exc_min", "exc_max", "ymin", "ymax"))
      }
    }else{ # if the regression didn't work, set to NA all coeffs
      resultmean[, c(paste0(signalcolumn, c("_dampingtime", "_exccoeff", "_eqvalue"))) := NA]
      summary <- model
      estimated <- NULL
    }

    # Output the results for the function
    return(list(data = data, resultmean = resultmean, regression = summary, estimated = estimated))

  }

}
# doremi_generate_order1 ----------------------------------------------------
#' Generation of the first order differential equation solution
#'
#' \code{doremi_generate_order1} returns a vector containing the convolution of the Green equation and the excitation function.
#'
#' @param dampingtime Signal damping coefficient. It represents the characteristic response time of the solution of the differential equation.
#' It should be a signed integer (dampingtime>0), otherwise a warning will be generated.
#' @param inputvec Is a vector containing the values of the excitation signal.
#' @param inputtim Is a vector containing the time values of the excitation signal.
#' @keywords Green, excitation, first order differential equation
#' @return Returns a list containing two elements:
#' y is a vector containing the values calculated from the convolution of the Green function and the excitation vector.
#' t is a vector containing the time values
#' @examples
#' doremi_generate_order1(10,rep(c(0,1,0),c(20,30,50)),1:100)
#'
#' \dontrun{
#' doremi_generate_order1(20,matrix(1,2,3),1:100)
#' }
#'@export
#'@importFrom stats convolve
doremi_generate_order1 <- function(dampingtime,
                                   inputvec,
                                   inputtim){
  #Error management
  #If inputvec is a scalar, the function warns the user that it should be a vector containing the values of the excitation signal
  if (length(inputvec) <= 1 | length(inputvec) != length(inputtim)) {
    warning("Both the excitation (inputvec) and is time values (inputtim) should be vectors and have the same length.\n")
  }
  #if inputvec is a character or a matrix, the function stops
  if (is.matrix(inputvec) | is.character(inputvec)) stop("inputvec should be a vector.")

  green <- exp(-1 / dampingtime * (inputtim - inputtim[1])) #Calculation of green function values.inputtim vector is forced to start in 0

  #as convolution doesn't take time into account, only number of points, a coefficient is needed to reduce amplitude
  convol <- (max(inputtim) - min(inputtim)) / (length(inputtim) - 1) * rev(convolve(rev(green), inputvec, type = "open"))
  convol <- convol[1:length(inputvec)]

  return(list(y = convol, t = inputtim))              #Returns the result of the convolution in the points corresponding to the excitation time vector.
}

# simulation_generate_order1 ----------------------------------------------

#' Simulation of various individual signals with intra and inter noise
#'
#' \code{simulation_generate_order1} This function generates signals with intra and inter noise for several individuals.

#' In order to do this, the function generates a pseudo-continuous signal per individual that is a solution to the first order differential equation
#' \deqn{\frac{dy(t)}{dt} - \gamma y(t) = E(t)}.
#' The analytical solution to this equation is a convolution between the Green function and the excitation (see function \code{\link{doremi_generate_order1}})
#' The function generates a pseudo-continuous signal to increase the precision with which the convolution is calculated. The number of points of this signal
#' are calculated with the inputs "precision" and "nf". From this expanded signal, the function samples points with a constant time step in order to have
#' the number of points indicated in nf. These operations are repeated as many times as the value set in the input "nindividuals". Once the signal is sampled,
#' intra and internoise with normal distributions are added.

#' @param nindividuals  Is a positive integer indicating the number of signals that will be simulated, each of them corresponding to an individual.
#' @param dampingtime Signal damping time. It is the characteristic response time of the solution to equation (1), and corresponds to the time needed to reach
#' 37\% of the maximum value when there is no excitation (or 63\% of the maximum value for a constant excitation). It should be a signed integer (dampingtime>0).
#' @inheritParams excitation_function
#' @param internoise Is the inter-individual noise added. The dampingtime accross individuals follows a normal distribution centered on the input parameter dampingtime
#' with a standard deviation of internoise*dampingtime
#' @param intranoise It is the noise added in each individual signal. It also follows a normal distribution with a standard deviation equal to this parameter times the maximum
#' amplitude of the convolution when using an excitation containing a single pulse (normalisation).

#' @keywords simulation, first order, differential equation
#' @return Returns a list of two objects:
#' rawdata- contains the function and time values evaluated in nf*precision points. This is to build a pseudo-continuous function with the parameters specified.
#' data- containts a sampling of raw data at equal time steps to keep only the number of points set in nf.
#' Each one of these elements contains the following columns:
#' id - individual identifier (starting by 1 and increasing by 1).
#' excitation - excitation signal generate through the excitation_function
#' timecol - time values
#' dampedsignalraw - signal with no noise (inter noise added for each individual)
#' dampedsignal - signal with intra noise added
#' @details Used for simulations in the context of the package.
#' @seealso \code{\link{doremi_generate_order1}} for calculation of the analytical solution to the differential equation
#' and \code{\link{excitation_function}} for excitation signal generation
#' @examples
#' simulation_generate_order1(nindividuals = 5,
#'                            dampingtime = 10,
#'                            amplitude = c(5,10),
#'                            nexc = 2,
#'                            duration = 20,
#'                            deltatf = 0.5,
#'                            tmax = 200,
#'                            minspacing = 0,
#'                            internoise = 0.2,
#'                            intranoise = 0.1)
#'@export
#'@importFrom data.table setDT
#'@importFrom data.table copy
#'@importFrom data.table :=
#'@importFrom data.table .N
#'@importFrom stats rnorm
simulation_generate_order1 <- function(nindividuals = 1,
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
  if (dampingtime<10) deltat=dampingtime/10
  else deltat=0.01

  if((deltatf < deltat) != 0){stop("Invalid deltatf, please modify.\n")}

  npoints <- tmax / deltat + 1

  # Generate simulation data for a given excitation and damping time
  # Generates id column (id being the individual number)
  # Creates a data table with id being the first column containing as many lines per individual as time points marked in npoints

  data <- setDT(list(id = unlist(lapply(c(1:nindividuals), function(x){rep(x, npoints)}))))

  #Add to that data table an "excitation" column, containing the values of the excitation signal.
  #Creates a new excitation signal for each individual
  data[, excitation := excitation_function(amplitude, nexc, duration, deltat, tmax, minspacing)$y, by = id]
  data[, timecol := excitation_function(amplitude, nexc, duration, deltat, tmax, minspacing)$t, by = id]

  #Creates a damping time vector by taking dampingtime input value and adding the internoise in a normal distribution
  #Contains as many elements as individuals
  dampingtimevec <- dampingtime + rnorm(nindividuals, mean = 0, sd = internoise * dampingtime)

  #If any value of the damping time vector is negative, the original value is used instead at that position of the vector
  if (any(dampingtimevec < 0)){
    dampingtimevec[dampingtimevec < 0] <- dampingtime
    warning("Some values for dampingtime where negative when adding internoise. The original dampingtime value has been used instead for those cases")
  }

  #Creates the signals for each individual taking the damping time for that individual from dampingtimevec
  #dampedsignalraw is the signal WITHOUT NOISE
  data[, dampedsignalraw := doremi_generate_order1 (dampingtimevec[.GRP], excitation, timecol)$y, by = id ]


  #Creates the signal for each individual with intra noise
  #In order to avoid adding an increased intra noise that will "grow" with each new pulse of the excitation signal, first
  #a "amplitudenorm" function is calculated in which the excitation has only one pulse that starts at the begining of the time sequence
  #The intranoise is added regarding the maximum value of that function
  #As the excitation function can receive an amplitude vector and a duration vector, in order to normalize, it will
  #be necessary to pick the maximum of the amplitude and the excitation
  if (length(amplitude) > 1){amp <- max(amplitude)}
  else {amp <- amplitude}
  if (length(duration) > 1){dur <- max(duration) / deltat + 1} #As done in the excitation function
  else{dur <- duration / deltat + 1}

  data[, amplitudenorm := doremi_generate_order1(dampingtimevec[.GRP], rep(c(amp, 0), c(dur, (length(excitation)-dur))), timecol)$y, by = id ]
  data[, dampedsignal := dampedsignalraw + rnorm(.N, mean = 0, sd = intranoise * max(abs(amplitudenorm))), by = id ]

  #Returning rawdata
  rawdata <- copy(data)

  #Selecting equally spaced data with deltatf
  #definition of time step for equally spaced values.
  data<-data[timecol %in% seq(0,tmax,deltatf), .SD,by = id] #Keeping only values from data table that are in the time intervals chosen

  return(list(rawdata = rawdata[, !"amplitudenorm", with = FALSE], data = data[, !c("amplitudenorm"), with = FALSE]))
}
