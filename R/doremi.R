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

# calculate.glla ----------------------------------------------------------
#' Calculation of derivatives using the GLLA method
#'
#' \code{calculate.glla} estimates the derivatives of a variable using the Generalized Local Linear Approximation (GLLA) method
#' described in \href{https://doi.org/10.4324/9780203864746}{Boker et al.(2010)}.
#' This method allows to estimate the derivatives over a number of measurement points called the embedding number.
#' @param signal is the input vector containing the data from which the derivatives are estimated.
#' @param time is a vector containing the time values corresponding to the signal. Arguments signal and time must have the same length.
#' @param embedding is an integer indicating the embedding dimension, that is the number of points to consider for derivative calculation.
#' Embedding must be at least #' 2 for the calculation of the first derivative (first order models) and at least 3 for the calculation of
#' the second derivative (second order models).
#' @keywords derivative, embedding dimension, rollmean
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
#' result <- calculate.glla(signal = signal, time = time, embedding = 5)
#'
#'@export
#'@importFrom zoo rollmean
calculate.glla <-  function(signal,
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
      deltat <- diff(time)[1]

      L1 <- rep(1,embedding)
      L2 <- c(1:embedding)-mean(c(1:embedding))
      L3 <- (L2^2)/2
      L <- cbind(L1,L2,L3)
      W <- L%*%solve(t(L)%*%L)

      Xembed <- embed(signal, embedding)
      Xembed <- Xembed[, ncol(Xembed):1]

      derivative <- Xembed %*% W
      derivative[,2] <- derivative[,2]/deltat
      derivative[,3] <- derivative[,3]/deltat^2
      derivative <- rbind(derivative, matrix(data = NA, ncol = 3, nrow = embedding - 1))

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
# calculate.fda ----------------------------------------------------------
#' Calculation of derivatives using the FDA method (splines)
#'
#' \code{calculate.fda} estimates the derivatives of a variable using the Functional Data Analysis (FDA)
#' method described in several sources, such as in \href{ISBN 978-0-387-98185-7}{Ramsay et al. (2009)}
#' and  \href{https://doi.org/10.1080/00273171.2015.1123138}{Chow et al. (2016)}.
#' This method estimates a spline function that fits all the data points and then derivates this function to estimate derivatives at those points.
#' In order for the derivatives to exist, the function must be smooth. A roughness penalty function controlled by a smoothing parameter is then used
#' The estimations are done by using the R's "fda" package.
#' @param signal is a vector containing the data from which the derivative is estimated.
#' @param time is a vector containing the time values corresponding to the signal. Arguments signal and time must have the same length.
#' @keywords derivative, embed, rollmean, fda, spline
#' @return Returns a list containing three columns:
#'
#' dtime- contains the initial time values provided.
#'
#' dsignal- is a data.frame containing three columns and the same number of rows as the signal.
#' The first column is the signal data points, the second is the first derivative evaluated at those points,
#' and the third is the second derivative evaluated at those points.
#' @examples
#' #In the following example the derivatives for the function y(t) = t^2 are calculated.
#' #The expected results are:
#' #y'(t) = 2t and y''(t) = 2
#' time <- c(1:500)/100
#' signal <- time^2
#' result <- calculate.fda(signal = signal, time = time)
#'
#'@export
calculate.fda <-  function(signal,
                           time,
                           spar = NULL){
  f <- smooth.spline(time,signal,spar = spar)
  derivative <- cbind(predict(f)$y,predict(f,deriv = 1)$y,predict(f,deriv = 2)$y)
  returnobject <- list("dtime" = time,
                       "dsignal" = derivative,
                       "spar" = spar)
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
#' @param duration is a vector of values greater or equal to 0 indicating the duration of each pulse in time units. It should have as many elements as the number of pulses (nexc). If
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
  if (any(duration < 0)){
    stop("Invalid input parameters. Pulse duration must be greater or equal to 0.\n")
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
# generate.1order ----------------------------------------------------
#' Generation of the first order differential equation solution with deSolve
#'
#' \code{generate.1order} returns a data frame containing the time (supplied as input) and a simulated signal generated as a solution to a first order
#' differential equation which coefficients are provided as inputs. The excitation is also provided as input and it can be null (then the solution
#' will be a decreasing exponential)
#' @param time Is a vector containing the time values corresponding to the excitation signal.
#' @param excitation Is a vector containing the values of the excitation signal.
#' @param y0 Signal initial value y(t=0)
#' @param tau Signal damping time. It represents the characteristic response time of the solution of the differential equation.
#' A negative value will produce divergence from equilibrium.
#' @param k Signal gain. It represents the proportionnality between the equilibrium value and the input maximum amplitude. It is thus relevant only
#' for differential equations including an excitation term.
#' @param yeq Signal equilibrium value. Value reached when the excitation term is 0 or constant.
#' @keywords first order differential equation, constant coefficients
#' @return Returns a list containing two elements:
#' \itemize{
#'   \item  y is a vector containing the values calculated with deSolve so that y is a solution to a first order differential equation with constant
#'   coefficients (provided as input).
#'   \item  t is a vector containing the corresponding time values
#' }
#' @examples
#' generate.1order(time = 0:49, excitation = c(rep(0,10),rep(1,40)))
#'@export
#'@importFrom deSolve ode
generate.1order <- function(time = 0:100,
                            excitation = NULL,
                            y0 = 0,
                            tau = 10,
                            k = 1,
                            yeq = 0){
  #Error management
  #If excitation is not supplied, then creation of an empty vector
  if(is.null(excitation)){excitation <- rep(0,length(time))}

  #If excitation is a scalar, the function warns the user that it should be a vector containing the values of the excitation signal
  if (length(excitation) <= 1 | length(excitation) != length(time)) {
    stop("Both the excitation (excitation) and its time values (time) should be vectors and have the same length.\n")
  }
  #if excitation is a character or a matrix, the function stops
  if (is.matrix(excitation) | is.character(excitation)) stop("Excitation should be a vector.")

  #Interpolating excitation function so that it can be evaluated in the time points required by deSolve. Rule 2 means constant interpolation.
  excf <- approxfun(time, excitation, rule = 2)

  #Initial values
  state <- c(y = y0)
  #Parameters
  parameters <- c(tau,k,yeq)

  #Model
  model1<-function(t, state, parameters)
  {   with(as.list(c(state, parameters)),
           { u <- excf(t)
           # return latent variable
           list(-(1/tau)*y + k/tau*u + yeq/tau)})

  }
  out <- as.data.frame(ode(y = state, times = time, func = model1, parms = parameters))
  names(out)<-c("t","y")
  out
}
# generate.2order ----------------------------------------------------
#' Generation of the second order differential equation solution with deSolve
#'
#' \code{generate.2order} returns a data frame containing the time (supplied as input) and a simulated signal generated as a solution to a second order
#' differential equation which constant coefficients are provided as inputs. The excitation is also provided as input and it can be null (then the solution
#' will be a damped linear oscillator)
#'
#' @param time is a vector containing the time values corresponding to the excitation signal.
#' @param excitation Is a vector containing the values of the excitation signal.
#' @param y0 is the initial condition for the variable (0, by default), it is a scalar.
#' @param v0 is the initial condition of the derivative (0, by default), it is a scalar.
#' @param xi is the damping factor. A negative value will produce divergence from equilibrium.
#' @param wn is the natural frequency which with the system would vibrate if there were no damping.
#' @param k is the gain. It represents the proportionnality between the equilibrium value and the input maximum amplitude. It is thus relevant only
#' for differential equations including an excitation term.
#' @param yeq is the signal equilibrium value. Value reached when the excitation term is 0 or constant.

#' @keywords second order differential equation, constant coefficients
#' @return Returns a list containing two elements:
#' \itemize{
#'   \item  y is a vector containing the values calculated with deSolve so that y is a solution to a second order differential equation with constant
#'   coefficients (provided as input).
#'   \item  t is a vector containing the corresponding time values
#' }
#' @examples
#' generate.2order(time=0:249,excitation=c(rep(0,10),rep(1,240)),wn=0.2)
#' generate.2order(y0=10)
#'@export
#'@importFrom deSolve ode
generate.2order <- function(time = 0:100,
                            excitation = NULL,
                            y0 = 0,
                            v0 = 0,
                            xi = 0.1,
                            period = 0.5,
                            k = 1,
                            yeq = 0){
  #Error management
  #If excitation is not supplied, then creation of an empty vector
  if(is.null(excitation)){excitation <- rep(0,length(time))}

  #If excitation is a scalar, the function warns the user that it should be a vector containing the values of the excitation signal
  if (length(excitation) <= 1 | length(excitation) != length(time)) {
    stop("Both the excitation (excitation) and its time values (time) should be vectors and have the same length.\n")
  }
  #if excitation is a character or a matrix, the function stops
  if (is.matrix(excitation) | is.character(excitation)) stop("Excitation should be a vector.")
  excf <- approxfun(time, excitation, rule = 2)

  #Initial values
  state <- c(y1 = y0, y2 = v0)
  parameters<-c(xi,period,k,yeq)

  #Model
  model2<-function(t, state, parameters)
  {   with(as.list(c(state, parameters)),
           { # equations of the state space model (two first oder coupled ODEs)
             u <- excf(t)
             wn <- 2*pi/period
             dy1 <- y2
             dy2 <- -wn^2*y1-2*xi*wn*y2 + k*wn^2*u + wn^2*yeq
             # return latent variables
             list(c(dy1,dy2))})
  }
  out <- as.data.frame(ode(y = state, times = time, func = model2, parms = parameters))
  names(out)<-c("t","y","dy")
  out
}
# generate.panel.1order ----------------------------------------------

#' Generation of first order differential equation solutions for several individuals with intra and inter noise
#'
#' \code{generate.panel.1order} Generation of first order differential equation solutions for several individuals with intra and inter noise.

#' The signals generated are solution to the first order differential equation:
#' \deqn{\frac{dy(t)}{dt} - \gamma y(t) = k*u(t) + yeq}
#' The analytical solution is generated with deSolve. From this
#' signal, the function samples points with a constant time step given by deltatf. These operations are repeated as many times as the value set in the input "nind". Once the signal is sampled,
#' intra-individual and inter-individual noise with normal distributions are added.

#' @inheritParams generate.1order
#' @param nind  number of individuals.
#' @param internoise Is the inter-individual noise added. The tau across individuals follows a normal distribution centered on the input parameter tau
#' with a standard deviation of internoise*tau, except if any damping time is negative (see Details section).
#' @param intranoise Is the noise added in each individual signal (dynamic noise). It also follows a normal distribution with a standard deviation equal to this parameter times the maximum
#' amplitude of the convolution when using an excitation containing a single pulse.

#' @keywords simulation, first order, differential equation
#' @return Returns a data frame with signal and time values starting at 0 and sampled at equal time steps deltatf in the time lapse tmax.
#' It contains the following columns:
#' \itemize{
#'    \item id - individual identifier (from 1 to nind).
#'    \item excitation - excitation signal generate through the generate.excitation
#'    \item time - time values
#'    \item signalraw - signal with no noise (inter noise added for each individual)
#'    \item dampedsignal - signal with intra noise added
#' }
#' @details Used for simulations in the context of the package.
#'
#' The function currently simulates only positive damping times corresponding to a regulated system. When the damping time is low
#' and the inter individual noise is high, some individuals' damping time could be negative. In that case, the damping time
#' distribution is truncated at 0.1*deltatf and values below are set to this limit. High values are symmetrically set at the upper percentile value
#' similar to a Winsorized mean. A warning provides the initial inter individual noise set as input argument and the inter individual
#' noise obtained after truncation.
#' @seealso \code{\link{generate.1order}} for calculation of the analytical solution to the differential equation
#' and \code{\link{generate.excitation}} for excitation signal generation
#' @examples
#' generate.panel.1order(time = generate.excitation(3, 6, 2, 1, 200, 2)$t,
#'                       excitation = generate.excitation(3, 6, 2, 1, 200, 2)$exc,
#'                       y0 = 0,
#'                       tau = 10,
#'                       k = 1,
#'                       yeq = 0,
#'                       nind = 5,
#'                       internoise = 0.2,
#'                       intranoise = 0.1)
#'@export
#'@importFrom data.table setDT
#'@importFrom data.table copy
#'@importFrom data.table :=
#'@importFrom data.table .N
#'@import stats
generate.panel.1order <- function(time,
                                  excitation = NULL,
                                  y0 = 0,
                                  tau = 10,
                                  k = 1,
                                  yeq = 0,
                                  nind = 1,
                                  internoise = 0,
                                  intranoise = 0){
  # Generating simulation data for a given excitation and parameters
  # Generating id column (id being the individual number)with as many lines per individual as npoints
  data <- data.table(id = rep(1:nind,each = length(time)))

  #Adding 3excitation and time to the data frame
  data[, time := time, by = id]
  data[, excitation := excitation, by = id]

  #Creating normal distributions for parameters
  tauvec <- rnorm(nind, mean = tau, sd = internoise * tau)
  y0vec <- rnorm(nind, mean = y0, sd = internoise * y0)
  kvec <- rnorm(nind, mean = k, sd = internoise * k)
  yeqvec <- rnorm(nind, mean = yeq, sd = internoise * yeq)

  #If any value of the damping time vector is negative, the original value is used instead at that position of the vector
  if (any(tauvec <= 0)){
    a <- 0.1 * deltatf #Threshold to truncate the tau distribution
    perc <- as.vector(prop.table(table(tauvec < a)))[2] #Calculation of the percentile of elements that are <0.1*deltatf (damping times of 0 are thus also excluded)
    b <- as.vector(quantile(tauvec, probs = 1 - perc)) #Calculation of the symmetrical threshold

    #Truncating the normal distribution to these two thresholds
    tauvec[tauvec < a] <- a
    tauvec[tauvec > b] <- b

    #Calculating new sd and internoise of the truncated distribution
    nsd <- sd(tauvec)
    ninternoise <- nsd / tau

    warning("Some values for tau where negative when adding internoise. Negative damping times imply signals
            that increase exponentially (diverge). This model generates signals that are self-regulated and thus,
            the normal distribution used has been truncated. The inter-individual noise added is of ", round(ninternoise,2), " instead of ",internoise,".\n")
  }
  #Addition of inter-noise
  #Creating the signals for each individual taking the parameters for that individual from the normal distribution vectors
  #signalraw is the signal WITHOUT NOISE
  data[, signalraw := generate.1order (time, excitation, y0vec[.GRP],tauvec[.GRP],kvec[.GRP],yeqvec[.GRP])$y, by = id ]

  #Addition of intra-noise
  data[, signal := signalraw + rnorm(.N, mean = 0, sd = intranoise * max(abs(signalraw))), by = id ]
  return(data)
}

# generate.panel.2order ----------------------------------------------

#' Generation of first order differential equation solutions for several individuals with intra and inter noise
#'
#' \code{generate.panel.2order} Generation of second order differential equation solutions for several individuals with intra and inter noise.

#' The signals are solution to the second order differential equation:
#' \deqn{\frac{d^2y}{dt} + 2\zeta\omega_{n}\frac{dy}{dt} + \omega_{n}^2 y = k*u(t)}
#' Where:
#' \itemize{
#'    \item{\eqn{\omega_{n} = \frac{2\pi}{period}} that is the system's natural frequency, corresponding to a period of period.
#'    The term \omega_{n}^2 represents thus the ratio between the attraction to the equilibrium and the inertia. If we considered the example
#'    of a mass attached to a spring, this term would represent the ratio of the spring constant and the object's mass.
#'    }
#'    \item{\zeta is the damping ratio. It represents the friction that damps the oscillation of the system (slows the rate of change of the variable).
#'    The term 2\zeta\omega_n thus represents the respective contribution of the inertia, the friction and the attraction to the equilibrium.
#'    The value of \zeta determines the shape of the system time response, which can be:
#'    \zeta<0	Unstable, oscillations of increasing magnitude
#'    \zeta=0	Undamped, oscillating
#'    0<\zeta<1	Underdamped or simply "damped". Most of the studies use this model, also referring to it as "Damped Linear Oscillator" (DLO).
#'    \zeta=1	Critically damped
#'    \zeta>1	Over-damped, no oscillations in the return to equilibrium
#'    }

#' The analytical solution is generated with deSolve. From this
#' signal, the function samples points with a constant time step given by deltatf. These operations are repeated as many times as the value set in the input "nind". Once the signal is sampled,
#' intra-individual and inter-individual noise with normal distributions are added.

#' @inheritParams generate.2order
#' #' @param nind  number of individuals.
#' @param internoise Is the inter-individual noise added. The damping ratio across individuals follows a normal distribution centered on the input parameter damping ratio
#' with a standard deviation of internoise*damping ratio, except if any damping time is negative (see Details section).
#' @param intranoise Is the noise added in each individual signal (dynamic noise). It also follows a normal distribution with a standard deviation equal to this parameter times the maximum
#' amplitude of the convolution when using an excitation containing a single pulse.

#' @keywords simulation, second order, differential equation
#' @return Returns a data frame with signal and time values for the time and excitation vectore provided.
#' It contains the following columns:
#' \itemize{
#'    \item id - individual identifier (from 1 to nind).
#'    \item excitation - excitation signal provided as input
#'    \item time - time values provided as input
#'    \item signalraw - signal with no noise (inter noise added for each individual)
#'    \item signal - signal with intra noise added
#' }
#' @details Used for simulations in the context of the package.
#'
#' The function currently simulates only positive damping factors corresponding to a self-regulated system. When the damping factor is low
#' and the inter individual noise is high, some individuals' damping factor could be negative. In that case, the damping factor
#' distribution is truncated at 0.1*min_deltat (calculated from the time vector) and values below are set to this limit.
#' High values are symmetrically set at the upper percentile value
#' similar to a Winsorized mean. A warning provides the initial inter individual noise set as input argument and the inter individual
#' noise obtained after truncation.
#' @seealso \code{\link{generate.2order}} for calculation of the analytical solution to the second order differential equation
#' and \code{\link{generate.excitation}} for excitation signal generation
#' @examples
#' generate.panel.2order(time = generate.excitation(3, 6, 2, 1, 200, 2)$t,
#'                       excitation = generate.excitation(3, 6, 2, 1, 200, 2)$exc,
#'                       y0 = 0,
#'                       v0 = 0,
#'                       xi = 0.1,
#'                       period = 0.5,
#'                       k = 1,
#'                       yeq = 0,
#'                       nind = 5,
#'                       internoise = 0.2,
#'                       intranoise = 0.1)
#'@export
#'@importFrom data.table data.table
#'@importFrom data.table copy
#'@importFrom data.table :=
#'@importFrom data.table .N
#'@import stats
generate.panel.2order <- function(time,
                                  excitation = NULL,
                                  y0 = 0,
                                  v0 = 0,
                                  xi = 0.1,
                                  period = 1,
                                  k = 1,
                                  yeq = 0,
                                  nind = 1,
                                  internoise = 0,
                                  intranoise = 0){
  # Generating simulation data for a given excitation and parameters
  # Generating id column (id being the individual number)with as many lines per individual as npoints
  data <- data.table(id = rep(1:nind,each = length(time)))

  #Adding excitation and time input vectors to the data frame
  data[, time := time, by = id]
  data[, excitation := excitation, by = id]


  #Creating normal distributions for parameters
  y0vec <- rnorm(nind, mean = y0, sd = internoise * y0)
  v0vec <- rnorm(nind, mean = v0, sd = internoise * v0)
  xivec <- rnorm(nind, mean = xi, sd = internoise * xi)
  periodvec <- rnorm(nind, mean = period, sd = internoise * period)
  kvec <- rnorm(nind, mean = k, sd = internoise * k)
  yeqvec <- rnorm(nind, mean = yeq, sd = internoise * yeq)

  #If any value of the damping time vector is negative, the original value is used instead at that position of the vector
  #this is because we are only interested in self-regulated signals, or oscillating but not divergent signals.
  if (any(xivec < 0)){
    a <- 0.1 * min(diff(time)) #Threshold to truncate the xi distribution
    perc <- as.vector(prop.table(table(xivec < a)))[2] #Calculation of the percentile of elements that are <0.1*deltatf (damping times of 0 are thus also excluded)
    b <- as.vector(quantile(xivec, probs = 1 - perc)) #Calculation of the symmetrical threshold

    #Truncating the normal distribution to these two thresholds
    xivec[xivec < a] <- a
    xivec[xivec > b] <- b

    #Calculating new sd and internoise of the truncated distribution
    nsd <- sd(xivec)
    ninternoise <- nsd / xi

    warning("Some values for xi where negative when adding internoise. Negative damping factors imply signals
            which oscillations increase exponentially (diverge). This function is intended to generate signals that are self-regulated and thus,
            the normal distribution used has been truncated. The inter-individual noise added is of ", round(ninternoise,2), " instead of ",internoise,".\n")
  }
  #Addition of inter-noise
  #Creating the signals for each individual taking the parameters for that individual from the normal distribution vectors
  #signalraw is the signal WITHOUT NOISE
  data[, signalraw := generate.2order(time, excitation, y0vec[.GRP],v0vec[.GRP],xivec[.GRP],periodvec[.GRP],kvec[.GRP],yeqvec[.GRP])$y, by = id ]

  #Addition of intra-noise
  #intranoise is the Signal to Noise ratio (SNR)
  #rnorm(signal) generates errors with mean 0 ans sd 1 then we need to calculate a scaling value
  #taken from https://stats.stackexchange.com/questions/31158/how-to-simulate-signal-noise-ratio
  if(intranoise!=0){
    data[, signal := signalraw + sqrt(var(signalraw)/(intranoise*var(rnorm(signalraw))))*rnorm(signalraw), by = id ]
  }else{
    data[, signal := signalraw, by = id]
  }
  return(data)
}

# analyze.1order --------------------------------------------------------------------
#' DOREMI first order analysis function
#'
#' \code{analyze.1order}  estimates the coefficients of a first order differential equation of the form:
#' \deqn{\frac{1}{\gamma} \dot{y}(t) = - y(t) + \epsilon u(t) + yeq}
#' using linear mixed-effect models.
#' Where y(t) is the individual's signal, \eqn{\dot{y}(t)} is the derivative and u(t) is the excitation.
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
#' @param derparam If dermethod "glla" or "gold" are chosen, it is the embedding number, a positive integer containing the number of points to be used for the calculation of the derivatives.
#' Its value by default is 2 as at least two points are needed for the calculation of the first derivative. If dermethod "fda" is chosen,
#' it is spar, the parameter related to the smoothing parameter lambda used in the penalization function of to estimate the derivatives via splines of \code{smooth.spline}
#' @param verbose Is a boolean that displays status messages of the function when set to 1.
#' @keywords analysis, first order, exponential
#' @return Returns a summary of the fixed components for the three coefficients: damping time, excitation coefficient and equilibrium value.
#' @details The analysis performs the following linear mixed-effects regression:
#' \deqn{y_{ij}'  \sim   b_{0} +b_{0j}+b_{1} y_{ij}+b_{2} E_{ij}+u_{1j} y_{ij}+u_{2j} E_{ij}+e_{ij}}
#' with i accounting for the time and j for the different individuals. \eqn{e_{ij}} are the residuals,
#' \eqn{y_{ij}'} is the derivative calculated on embedding points and
#' y and E are the signal and the excitation averaged on embedding points.
#' The coefficients estimated to characterize the signal are calculated as follows:
#' \itemize{
#'   \item Damping time, tau:  \eqn{\tau _{j} =  \frac{1}{ \gamma _{j} }}  with \eqn{\gamma _{j} =  b_{1} + u_{1j} }
#'   \item Gain, k: \eqn{\epsilon _{j} = \frac{b_{2} + u_{2j}}{\gamma _{j}}}. It is the proportionality between the excitation and the
#'   difference between the maximum value reached by the signal and its initial value.
#'   \item Equilibrium value, yeq: \eqn{yeq _{j} = \frac{b_{0} + b_{0j}}{\gamma _{j}}}. It is the stable value reached in the absence of excitation.
#' }
#' The estimation is performed using the function lmer if there are several individuals or lm if there is only one.
#' With the above estimated parameters, the estimated signal can be reconstructed for
#' each individual by using the generate.1order function on this package (based on deSolve's ode function).
#' The function returns five objects:
#' \enumerate{
#'  \item data- A data.frame including the input data, the intermediate calculations used to prepare the variables for
#'   the fit and the estimated trajectories for each individual.
#'
#'     signal_rollmean - calculation of the moving average of the signal over embedding points.
#'
#'     signal_derivate1 - calculation of the first derivative of the signal with the glla method in embedding points.
#'
#'     time_derivate - calculation of the moving average of the time vector over embedding points.
#'
#'     input_rollmean - calculation of the moving average of the excitation vector over embedding points.
#'
#'     estimated- the estimated signal calculated using deSolve's ode function with a first order model, the excitation provided as input and the damping time,
#'     excitation coefficient and equilibrium value obtained from the fit.
#'  \item resultid- A data.frame including for each individual, listed by id number, the damping time, the excitation coefficient and the
#'  equilibrium value (see variables presented in the Details section).
#'  \item resultmean- A data.frame including the fixed effects of the three coefficients mentioned above.
#'  \item regression- A list containing the summary of the linear mixed-effects regression.
#'
#'  As seen in the Description section, the print method by default prints only the resultmean element. Each one of the other objects
#'  can be accessed by indicating $ and their name after the result, for instance, for a DOREMI object called "result", it is possible
#'  to access the regression summary by typing result$regression.
#'  \item embedding - contains the embedding number used to generate the results (same as function input argument)
#'  \item spar - contains the smoothing parameter used for the estimation of the derivatives using splines (method "fda")
#' }
#' @seealso \code{\link{calculate.gold}\link{calculate.glla}\link{calculate.fda}} to compute the derivatives, for details on embedding/spar.
#' @examples
#' myresult <- analyze.1order(data = cardio,
#'                   id = "id",
#'                   input = "load",
#'                   time = "time",
#'                   signal = "hr",
#'                   dermethod ="calculate.gold",
#'                   derparam = 5)
#'@export
#'@import data.table
#'@importFrom lmerTest lmer
#'@importFrom lme4 lmerControl
#'@importFrom lme4 ranef
#'@importFrom stats embed
#'@importFrom stats lm
#'@importFrom zoo rollmean
#'@importFrom zoo na.locf
analyze.1order <- function(data,
                 id = NULL,
                 input = NULL,
                 time = NULL,
                 signal,
                 dermethod = "calculate.fda",
                 derparam = 2,
                 verbose = FALSE){

  intdata <- setDT(copy(data)) # Makes a copy of original data so that it can rename columns freely if needed. setDT converts it to data.table
  noinput <- FALSE #Flag that will allow to differentiate if there is an excitation term or not when doing the regression

# Error management
    errorcheck(intdata,signal)
  if (!is.null(id)){ #Several individuals
      errorcheck(intdata,id)
      nind <- length(unique(intdata[[id]]))
  }else{#Single individual. Create id column anyways (needed for the rest of the data processing)
      nind <-1
      intdata[, id:=1]
      id<- "id"
  }
  if (!is.null(input)){
    errorcheck(intdata, input)
    for (i in 1:length(input)){ #For loop to go through all the inputs
      if(all(intdata[, diff(get(input)), by = id]$V1 == 0) == TRUE){
        noinput <- TRUE
        warning("Excitation signal is constant. Adjustment will consist on exponential fit.\n")
      } #If input is constant
    }
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
  if(!dermethod %in% c("calculate.fda","calculate.glla","calculate.gold")){
    stop("Derivative method is not valid. Please introduce the name of the derivative calculation function: \"calculate.fda\",\"calculate.glla\" or \"calculate.gold\"")
  }

  #After verification, extracting only relevant columns
  intdata <- intdata[,.SD,.SDcols = c(id, input, time, signal)]

  #Sorting table
  setkeyv(intdata,c(id,time))

  #Find time duplicates and display error message if it is the case
  intdata[,timedup:=lapply(.SD,duplicated),.SDcols = time,by = id]
  if(any(intdata$timedup)){stop("Input data.table contains duplicated time points.\n")}
    else  intdata[, timedup := NULL]

  #Verifying column names repeated in data table.
  if(any(duplicated(colnames(intdata)))){stop("Input datatable contains duplicated column names.\n")}

  #Suppress NA elements if there is an NA in time or in signal
  #This suppresses the entire row of the intdata table
  intdata <- intdata[!is.na(time) & !is.na(signal)]

  #If in the intdata left, there are some NA values in the excitation vector, set them to zero.
  #This is to avoid losing intdata from signal and time vectors, as sometimes when people fill in the signal intdata they put NA to mean no excitation, thus 0.
  na.to.0 <- function(x){x[is.na(x)] <- 0; x}
  intdata[, (input) := lapply(.SD, na.to.0), .SDcols = input]

  #Saving original names and then renaming data table to internal names
  originalnames <- c(id, time, signal, input)
  doremiexc <- paste0("input", seq(input)) # Doremi excitation vector ("input1","input2","input3"...)
  doreminames <- c("id_tmp", "time", "signal", doremiexc)
  setnames(intdata, originalnames, doreminames)
  #Rename id column to "id_tmp" and create an extra column "id" that is numeric
  #Regression works better with numeric id (comes from definition of lmer)
  intdata[, id := rep(1:length(unique(id_tmp)), intdata[, .N, by = id_tmp]$N)]


  #Calculation of the signal rollmean and first derivative of the signal column
  derivate <-get(dermethod)
  intdata[, signal_rollmean := derivate(signal, time, derparam)$dsignal[, 1], by = id]
  intdata[, signal_derivate1 := derivate(signal, time, derparam)$dsignal[, 2], by = id]
  intdata[, time_derivate := derivate(signal, time, derparam)$dtime, by = id]

  #Calculation of the roll mean of the excitation columns if there is at least one input column
  if (!noinput){
    #myfun <- function(x){x[] <- c(rollmean(x, (derparam)), rep(NA, derparam - 1)); x}
    #intdata[, (paste0(doremiexc,"_rollmean")) := lapply(.SD, myfun), .SDcols = doremiexc, by = id]
    intdata[, (paste0(doremiexc,"_rollmean")) := lapply(.SD, function(x){x[] <- derivate(x,time,derparam)$dsignal[,1];x}), .SDcols = doremiexc, by = id]
  }

  #Linear mixed-effect regression MULTIPLE INDIVIDUALS
  if(nind>1){
    if (noinput){ # if there is no excitation signal
      model <- tryCatch({lmer(signal_derivate1 ~ signal_rollmean + (1 + signal_rollmean |id),
                              data = intdata, REML = TRUE, control = lmerControl(calc.derivs = FALSE, optimizer = "nloptwrap"))}, error = function(e) e)
      if (verbose){print("Status: Unknown excitation. Linear mixed-effect model calculated.")}
    }else{ # if there is one OR SEVERAL excitation signals
      model <- tryCatch({lmer(paste0("signal_derivate1 ~ signal_rollmean + (1 +", paste(doremiexc, "rollmean ", collapse = "+", sep = "_"),
                              " + signal_rollmean |id) + ", paste(doremiexc, "rollmean ", collapse = "+",sep = "_")),
                              data = intdata, REML = TRUE, control = lmerControl(calc.derivs = FALSE, optimizer = "nloptwrap"))}, error = function(e) e)
      if (verbose){print("Status: One or several excitations. Linear mixed-effect model calculated.")}
    }
  }else{ #SINGLE individual
    if(noinput){ # if there is no excitation signal
      model <- tryCatch({lm(signal_derivate1 ~ signal_rollmean, data = intdata)}, error = function(e) e)
      if (verbose){print("Status: Unknown excitation. Linear regression calculated")}

    }
    else{ # if there is one or several excitation signals
      model <- tryCatch({lm(paste0("signal_derivate1 ~ signal_rollmean + ", paste(doremiexc, "rollmean ", collapse = "+", sep = "_")),
                            data = intdata)}, error = function(e) e)
      if (verbose){print("Status: One or several excitations. Linear regression calculated")}
    }
  }
    if (!inherits(model,"error")){ # if the regression worked
      if (verbose){print("Status: Linear mixed-effect model had no errors.")}
      summary <- summary(model) # Summary of the regression
      if(nind > 1){random <- ranef(model)} # Variation of the estimated coefficients over the individuals. Only for data with several individuals
      if(nind > 1){regression <- list(summary, random)}else{regression <- summary} # list to output both results: summary, and the table from ranef

      #Create tables with results
      #The first table contains the input data

      #The second table contains the mean values for gamma and thao for all the individuals (single line)
      #Generate mean results with convergence criterions
      resultmean <- setDT(list(id = "All"))

      # calculate the damping time for all signal columns: -1/damping_coeff
      resultmean[, tau := -1L/summary$coefficients["signal_rollmean", "Estimate"]]

      # Extract the intercept coeff (equilibrium value)
      resultmean[, yeq := summary$coefficients["(Intercept)","Estimate"] * resultmean[, tau]]

      if(nind > 1){
        #The third table contains the results for gamma and thao for each individual (one line per individual)
        resultid <- setDT(list(id = unique(intdata$id), id_tmp = unique(intdata$id_tmp)))
        setkey(resultid, id) #sorts the data table by id

        # Damping time in resultid
        resultid[, tau := -1L/(summary$coefficients["signal_rollmean","Estimate"] + random$id[.GRP,"signal_rollmean"]), by = id]

        # Extract the intercept (equilibrium value) calculated for each individual (present in random, regression table)
        # Offset in resultid
        resultid[, yeq := (summary$coefficients["(Intercept)", "Estimate"] + random$id[.GRP, "(Intercept)"]) * resultid[.GRP, tau], by = id]

        # Generation of the estimated signal for all id using analyze.1order generate FOR SEVERAL INDIVIDUALS (will be added to the $data object)
        # Write a warning if any of the damping times calculated was negative
        if (any(is.na(resultid[, tau])) | any(resultid[, tau] < 0)){
          warning("Some of the damping times calculated were negative and thus, the estimated signal was not generated for these.
                  Damping times can be negative for some individuals for the following reasons: 1. The signal of
                  the individual doesn't go back to equilibrium. 2.The linear mixed-effects model estimating the random
                  effect showed some error messages/warnings. 3.Model misspecification.\n")
        }
      }else{resultid <- NULL} #Single individual will not have resultid table

      if (noinput){ #There is no excitation signal as input: exponential fit
        if (verbose){print("Status: Unknown excitation. Calculation of estimated signal through exponential fit.")}
        #Exponential model is assumed when no input is provided and thus a new fit is necessary BY INDIVIDUAL
        #log(y-B) = gamma*t + log(A) --> LINEAR EQUATION
        #B is known, it is the intercept*tau, from the model fit with the derivative inside the analysis function
        #A is unknown. It can be found by fitting log(y-B)~ t
        #We'll call the coefficients resulting from this new fit Ap and Bp: Ap=gamma, Bp=log(A)
        if(nind > 1){ # Several individuals
          intdata[, expfit_A := lm(log(abs(signal - resultid[.GRP, yeq])) ~ time)$coefficients[1], by = id]
          intdata[, signal_estimated :=
                    if(!is.na(resultid[.GRP, tau]) && resultid[.GRP, tau] > 0){
                      # if there is a damping time that has been calculated and if it is greater than 0 (decreasing exponential)
                      #Then it is assumed that the signal follows a decreasing exponential: y = A* exp(gamma*t)+B
                      #expmodel comes from fitting log(y-B)~t. Calculated above
                      exp(expfit_A) * exp(-1L / resultid[.GRP, tau] * time) + resultid[.GRP, yeq]
                    }else{NaN}, by = id]
        }else{ #Single individual
          intdata[, expfit_A := lm(log(abs(signal - resultmean[, yeq])) ~ time)$coefficients[1]]
          intdata[, signal_estimated :=
                    if(!is.na(resultmean[, tau]) && resultmean[, tau] > 0){
                      # if there is a damping time that has been calculated and if it is greater than 0 (decreasing exponential)
                      #Then it is assumed that the signal follows a decreasing exponential: y = A* exp(gamma*t)+B
                      #expmodel comes from fitting log(y-B)~t. Calculated above
                      exp(expfit_A) * exp(-1L / resultmean[, tau] * time) + resultmean[, yeq]
                    }else{NaN}]
        }
        #Removing temporary column "expfit_A" from table "intdata"
        intdata <- intdata[, c("expfit_A") := NULL]

      }else{
        if (verbose){print("Status: One or several excitation terms. Calculation of estimated signal with deSolve")}
        # Extract the excitation coeff for each excitation
        intdata[, totalexc := 0]
        intdata[, totalexcroll := 0]
        for (i in 1:length(input)){ #For loop to go through all the inputs
          resultmean[, paste0(doremiexc[i],"_k") := summary$coefficients[paste0(doremiexc[i], "_rollmean"), "Estimate"] * resultmean[, tau]]

          #If variation of the excitation coefficient across individuals needed:
          #And for each individual: the mean coeff (sumary$coeff) + the variation per Individual (in random)
          #Excitation coefficient in resultid
          if(nind > 1){  #Several individuals
            resultid[, paste0(doremiexc[i],"_k") := (summary$coefficients[paste0(doremiexc[i], "_rollmean"), "Estimate"] +
                                                           random$id[.GRP,paste0(doremiexc[i], "_rollmean")]) * resultid[.GRP, tau], by = id]
            intdata[, totalexc := totalexc +  (summary$coefficients[paste0(doremiexc[i], "_rollmean"), "Estimate"] +
                                                 random$id[.GRP,paste0(doremiexc[i], "_rollmean")]) * resultid[.GRP, tau]* get(doremiexc[i]), by = id] #total excitation is the sum of all the excitations
            intdata[, totalexcroll := totalexcroll +  (summary$coefficients[paste0(doremiexc[i], "_rollmean"), "Estimate"] +
                                                 random$id[.GRP,paste0(doremiexc[i], "_rollmean")]) * resultid[.GRP, tau]* get(paste0(doremiexc[i],"_rollmean")), by = id]

            #ponderated with their corresponding gains
          }else{ #Single individual
            intdata[, totalexc := totalexc + summary$coefficients[paste0(doremiexc[i], "_rollmean"), "Estimate"] * resultmean[, tau] * get(doremiexc[i])]
            intdata[, totalexcroll := totalexcroll + summary$coefficients[paste0(doremiexc[i], "_rollmean"), "Estimate"] * resultmean[, tau] * get(paste0(doremiexc[i],"_rollmean"))]

          }
        }
        #The estimated signal is calculated by calling ode function in deSolve (through function "generate.1order"). As we will have a decomposition
        #of k for each excitation, the excitation considered is already the total excitation with the total gain (to avoid calculating both separately, this is why
        #k=1, total gain is already included in totalexc)
        #Assuming initial value is equilibrium value
        generate_solution<- function(time,time_derivate,totalexc,totalexcroll,signal_rollmean,tau,yeq){
          #Building time vector that had the first point of rollmean time and rollmean signal
          timecomp<-c(time[time<time_derivate[1]],time_derivate[1],time[time>time_derivate[1]])
          exccomp<-c(totalexc[time<time_derivate[1]],totalexcroll[1],totalexc[time>time_derivate[1]])
          #Position of the rollmean time in the original time vector
          i <- length(time[time<time_derivate[1]])+1
          droite <- generate.1order(time = timecomp[i:length(timecomp)], #includes original time and the first element of the rollmean time
                                    excitation = exccomp[i:length(exccomp)],
                                    y0 = signal_rollmean[1],
                                    tau = tau,
                                    k = 1,
                                    yeq = yeq)
          gauche <- generate.1order(time = timecomp[i:1],
                                    excitation = exccomp[i:1],
                                    y0 = signal_rollmean[1],
                                    tau = tau,
                                    k = 1,
                                    yeq = yeq)
          #putting the two solutions together
          sol<-c(gauche$y[i:2],droite$y)
          sol
        }
        if(nind > 1){
          intdata[, signal_estimated := generate_solution(time,time_derivate,totalexc,totalexcroll,signal_rollmean,resultid[.GRP, tau],resultid[.GRP, yeq]),by =id]
        }else{
          intdata[, signal_estimated := generate_solution(time,time_derivate,totalexc,totalexcroll,signal_rollmean,resultmean[, tau],resultmean[, yeq])]
        }
      }

    }else{ # if the regression didn't work, a warning will be generated and tables will be set to NULL
      if (verbose){print("Status: Linear mixed-effect model produced errors.")}
      warning("Linear mixed-effect regression produced an error. Verify the regression object of the result.\n")
      resultid <- NULL
      resultmean <- NULL
      regression <- model
    }

    #Renaming columns in $data, $resultid, $resultmean objects to original names
    intdata[, id := NULL]
    if(!is.null(resultid)){resultid[, id := NULL]}

    #Replacing names in $data
    intdatanames <- names(intdata)
    intdatanamesnew <- names(intdata)

    rmeannames <- names(resultmean)
    rmeannamesnew <- names(resultmean)

    ridnames <- names(resultid)
    ridnamesnew <- names(resultid)

    for(idx in seq(doreminames)){
      intdatanamesnew <-gsub(paste0("^",doreminames[idx],"(?![0-9])"), originalnames[idx], intdatanamesnew, perl = T)
      #Replacing names in $resultmean
      if(!is.null(resultid)){
        ridnamesnew <-gsub(paste0("^",doreminames[idx],"(?![0-9])"), originalnames[idx], ridnamesnew, perl = T)
      }
      if(!is.null(resultmean)){
        rmeannamesnew <-gsub(paste0("^",doreminames[idx],"(?![0-9])"), originalnames[idx], rmeannamesnew, perl = T)
      }
    }
    setnames(intdata, intdatanames, intdatanamesnew)
    if(!is.null(resultid)){setnames(resultid, ridnames, ridnamesnew)}
    if(!is.null(resultmean)){setnames(resultmean, rmeannames, rmeannamesnew)}

  #Output the results of the function
  #Excitation string
  #If there is no excitation term
  if (noinput){str_exc <- 0
  }else{str_exc <- input} # If there is one OR SEVERAL excitation columns
  res = list(data = intdata, resultid = resultid, resultmean = resultmean, regression = regression, dermethod = dermethod, derparam = derparam, str_time = time, str_exc = str_exc, str_signal = signal, str_id = id)
  class(res)= "doremi" # Class definition
  return(res)
}

# analyze.2order --------------------------------------------------------------------
#' DOREMI second order analysis function
#'
#' \code{analyze.2order}  estimates the coefficients of a second order differential equation of the form:
#' \deqn{\frac{d^2y}{dt} + 2\zeta\omega_{n}\frac{dy}{dt} + \omega_{n}^2 (y - y_{eq}) = k*u(t) }
#' Where y(t) is the individual's signal, \eqn{\dot{y}(t)} is the derivative and u(t) is the excitation.
#' The function estimates the coefficients \zeta, \omega_{n}, k and y_{eq} using a two step procedure in
#' which the user can choose the derivative estimation method (through the parameter dermethod) and then the
#' coefficients are then estimated via linear mixed-effect models.
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
#' @param embedding Is a positive integer containing the number of points to be used for the calculation of the derivatives if dermethod "glla" or "gold" are chosen.
#' (It will be ignored if method "fda" is chosen) Its value by default is 2 as at least two points are needed for the calculation of the first derivative.
#' @param spar Related to the smoothing parameter lambda used in the penalization function of to estimate the derivatives via splines. If dermethod "glla" or "gold"
#' are chosen, this parameter is ignored. For more details, see the documentation
#' of \code{smooth.spline}
#' @param verbose Is a boolean that displays status messages of the function when set to 1.
#' @keywords analysis, second order

#' @return Returns a summary of the fixed components for the coefficients: damping factor, natural frenquency, gain and equilibrium value.
#' The estimation is performed using the function lmer if there are several individuals or lm if there is only one.
#' With the above estimated parameters, the estimated signal is reconstructed for
#' each individual by using the generate.2order function on this package (based on deSolve's ode function).
#' The function returns five objects:
#' \enumerate{
#'  \item data- A data.frame including the input data, the intermediate calculations used to prepare the variables for
#'   the fit and the estimated trajectories for each individual.
#'
#'     signal_rollmean - calculation of the moving average of the signal over embedding points.
#'
#'     signal_derivate1 - calculation of the first derivative of the signal with the chosen method in embedding points/with smoothing parameter spar
#'
#'     signal_derivate2 - calculation of the second derivative of the signal with the chosen method in embedding points/with smoothing parameter spar
#'
#'     time_derivate - calculation of the moving average of the time vector over embedding points.
#'
#'     input_rollmean - calculation of the moving average of the excitation vector over embedding points.
#'
#'     estimated- the estimated signal calculated using deSolve's ode function with a second order model, the excitation provided as input and the
#'     coefficients obtained from the fit.
#'  \item resultid- A data.frame including for each individual, listed by id number, the coefficients calculated (thus fixed + random component)
#'  \item resultmean- A data.frame including the fixed effects of the coefficients mentioned above.
#'  \item regression- A list containing the summary of the linear mixed-effects regression.
#'
#'  As seen in the Description section, the print method by default prints only the resultmean element. Each one of the other objects
#'  can be accessed by indicating $ and their name after the result, for instance, for a DOREMI object called "result", it is possible
#'  to access the regression summary by typing result$regression.
#'  \item embedding - contains the embedding number used to generate the results (same as function input argument). Will only appear in the output if the
#'  derivative estimation method chosen is "glla" or "gold".
#'  \item spar -
#' }
#' @seealso \code{\link{calculate.gold}\link{calculate.glla}\link{calculate.fda}} to compute the derivatives, for details on embedding.
#' @examples
#' myresult <- analyze.2order(data = cardio,
#'                   id = "id",
#'                   input = "load",
#'                   time = "time",
#'                   signal = "hr",
#'                   dermethod = "calculate.gold",
#'                   derparam = 5)
#'@export
#'@import data.table
#'@importFrom lmerTest lmer
#'@importFrom lme4 lmerControl
#'@importFrom lme4 ranef
#'@importFrom stats embed
#'@importFrom stats lm
#'@importFrom zoo rollmean
#'@importFrom zoo na.locf
analyze.2order <- function(data,
                           id = NULL,
                           input = NULL,
                           time = NULL,
                           signal,
                           dermethod = "calculate.gold",
                           derparam = 3,
                           verbose = FALSE){

  intdata <- setDT(copy(data)) # Makes a copy of original data so that it can rename columns freely if needed. setDT converts it to data.table
  noinput <- FALSE #Flag that will allow to differentiate if there is an excitation term or not when doing the regression

  # Error management
  errorcheck(intdata,signal)
  if (!is.null(id)){ #Several individuals
    errorcheck(intdata,id)
    nind <- length(unique(intdata[[id]]))
  }else{#Single individual. Create id column anyways (needed for the rest of the data processing)
    nind <-1
    intdata[, id:=1]
    id<- "id"
  }
  if (!is.null(input)){
    errorcheck(intdata, input)
    for (i in 1:length(input)){ #For loop to go through all the inputs
      if(all(intdata[, diff(get(input)), by = id]$V1 == 0) == TRUE){
        noinput <- TRUE
        warning("Excitation signal is constant. Adjustment will consist on exponential fit.\n")
      } #If input is constant
    }
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
  if(!dermethod %in% c("calculate.fda","calculate.glla","calculate.gold")){
    stop("Derivative method is not valid. Please introduce the name of the derivative calculation function: \"calculate.fda\",\"calculate.glla\" or \"calculate.gold\"")
  }

  #After verification, extracting only relevant columns
  intdata <- intdata[,.SD,.SDcols = c(id, input, time, signal)]

  #Sorting table
  setkeyv(intdata,c(id,time))

  #Find time duplicates and display error message if it is the case
  intdata[,timedup:=lapply(.SD,duplicated),.SDcols = time,by = id]
  if(any(intdata$timedup)){stop("Input data.table contains duplicated time points.\n")}
  else  intdata[, timedup := NULL]

  #Verifying column names repeated in data table.
  if(any(duplicated(colnames(intdata)))){stop("Input datatable contains duplicated column names.\n")}

  #Suppress NA elements if there is an NA in time or in signal
  #This suppresses the entire row of the intdata table
  intdata <- intdata[!is.na(time) & !is.na(signal)]

  #If in the intdata left, there are some NA values in the excitation vector, set them to zero.
  #This is to avoid losing intdata from signal and time vectors, as sometimes when people fill in the signal intdata they put NA to mean no excitation, thus 0.
  na.to.0 <- function(x){x[is.na(x)] <- 0; x}
  intdata[, (input) := lapply(.SD, na.to.0), .SDcols = input]

  #Saving original names and then renaming data table to internal names
  originalnames <- c(id, time, signal, input)
  doremiexc <- paste0("input", seq(input)) # Doremi excitation vector ("input1","input2","input3"...)
  doreminames <- c("id_tmp", "time", "signal", doremiexc)
  setnames(intdata, originalnames, doreminames)

  #Rename id column to "id_tmp" and create an extra column "id" that is numeric
  #Regression works better with numeric id (comes from definition of lmer)
  intdata[, id := rep(1:length(unique(id_tmp)), intdata[, .N, by = id_tmp]$N)]


  #Calculation of the signal rollmean and first derivative of the signal column
  derivate<-get(dermethod)
  intdata[, signal_rollmean := derivate(signal, time, derparam)$dsignal[, 1], by = id]
  intdata[, signal_derivate1 := derivate(signal, time, derparam)$dsignal[, 2], by = id]
  intdata[, signal_derivate2 := derivate(signal, time, derparam)$dsignal[, 3], by = id]
  intdata[, time_derivate := derivate(signal, time, derparam)$dtime, by = id]

  #Calculation of the roll mean of the excitation columns if there is at least one input column
  if (!noinput){
    # myfun <- function(x){x[] <- c(rollmean(x, (derparam)), rep(NA, derparam - 1)); x}
    # intdata[, (paste0(doremiexc,"_rollmean")) := lapply(.SD, myfun), .SDcols = doremiexc, by = id]
    intdata[, (paste0(doremiexc,"_rollmean")) := lapply(.SD, function(x){x[] <- derivate(x,time,derparam)$dsignal[,1];x}), .SDcols = doremiexc, by = id]
  }

  #Linear mixed-effect regression MULTIPLE INDIVIDUALS
  if(nind>1){
    if (noinput){ # if there is no excitation signal
      model <- tryCatch({lmer(signal_derivate2 ~ signal_derivate1 + signal_rollmean + (1 + signal_rollmean + signal_derivate1 |id),
                              data = intdata, REML = TRUE, control = lmerControl(calc.derivs = FALSE, optimizer = "nloptwrap"))}, error = function(e) e)
      if (verbose){print("Status: Unknown excitation. Linear mixed-effect model calculated.")}
    }else{ # if there is one OR SEVERAL excitation signals
      model <- tryCatch({lmer(paste0("signal_derivate2 ~ signal_derivate1 + signal_rollmean + (1 +", paste(doremiexc, "rollmean ", collapse = "+", sep = "_"),
                                     " + signal_rollmean + signal_derivate1 |id) + ", paste(doremiexc, "rollmean ", collapse = "+",sep = "_")),
                                     data = intdata, REML = TRUE, control = lmerControl(calc.derivs = FALSE, optimizer = "nloptwrap"))}, error = function(e) e)
      if (verbose){print("Status: One or several excitations. Linear mixed-effect model calculated.")}
    }
  }else{ #SINGLE individual
    if(noinput){ # if there is no excitation signal
      model <- tryCatch({lm(signal_derivate2 ~ signal_derivate1 + signal_rollmean, data = intdata)}, error = function(e) e)
      if (verbose){print("Status: Unknown excitation. Linear regression calculated")}

    }
    else{ # if there is one or several excitation signals
      model <- tryCatch({lm(paste0("signal_derivate2 ~ signal_derivate1 + signal_rollmean + ", paste(doremiexc, "rollmean ", collapse = "+", sep = "_")),
                            data = intdata)}, error = function(e) e)
      if (verbose){print("Status: One or several excitations. Linear regression calculated")}
    }
  }
  if (!inherits(model,"error")){ # if the regression worked
    if (verbose){print("Status: Linear mixed-effect model had no errors.")}
    summary <- summary(model) # Summary of the regression
    if(nind > 1){random <- ranef(model)} # Variation of the estimated coefficients over the individuals. Only for data with several individuals
    if(nind > 1){regression <- list(summary, random)}else{regression <- summary} # list to output both results: summary, and the table from ranef

    #Create tables with results
    #The first table contains the input data
    #The second table contains the mean values for gamma and thao for all the individuals (single line)
    #Generate mean results with convergence criterions
    resultmean <- data.table(omega2 = -1*summary$coefficients["signal_rollmean","Estimate"],
                             omega2_std = summary$coefficients["signal_rollmean","Std. Error"],
                             esp2omega = -1*summary$coefficients["signal_derivate1","Estimate"],
                             esp2omega_std = summary$coefficients["signal_derivate1","Std. Error"],
                             yeqomega2 = summary$coefficients["(Intercept)","Estimate"],
                             yeqomega2_std = summary$coefficients["(Intercept)","Std. Error"])
    resultmean[,period := ifelse(omega2>0,2*pi/sqrt(omega2),NaN)]
    resultmean[, wn := ifelse(omega2>0,sqrt(omega2),NaN)]
    resultmean[, xi := ifelse(omega2>0,esp2omega/(2*sqrt(omega2)),NaN)]
    resultmean[, yeq := ifelse(omega2>0,yeqomega2/(omega2),NaN)]

    if(nind > 1){
      #The third table contains the results for the coefficients for each individual (one line per individual)
      resultid <- setDT(list(id = unique(intdata$id), id_tmp = unique(intdata$id_tmp)))
      setkey(resultid, id) #sorts the data table by id

      resultid[, omega2 := -(summary$coefficients["signal_rollmean","Estimate"] + random$id[.GRP,"signal_rollmean"]), by = id]
      resultid[, esp2omega := -(summary$coefficients["signal_derivate1","Estimate"] + random$id[.GRP,"signal_derivate1"]), by = id]
      resultid[, yeqomega2 := summary$coefficients["(Intercept)", "Estimate"] + random$id[.GRP, "(Intercept)"], by = id]

      resultid[,period := ifelse(omega2>0,2*pi/sqrt(omega2),NaN)]
      resultid[, wn := ifelse(omega2>0,sqrt(omega2),NaN)]
      resultid[, xi := ifelse(omega2>0,esp2omega/(2*sqrt(omega2)),NaN)]
      resultid[, yeq := ifelse(omega2>0,yeqomega2/(omega2),NaN)]

      # Generation of the estimated signal for all id using analyze.2order generate FOR SEVERAL INDIVIDUALS (will be added to the $data object)
      # Write a warning if any of the damping times calculated was negative
      if (any(is.na(resultid[, xi])) | any(resultid[, xi] < 0)){
        warning("Some of the damping factors calculated were negative and thus, the estimated signal was not generated for these.
                  Damping factors can be negative for some individuals for the following reasons: 1. The signal of
                  the individual doesn't go back to equilibrium. 2.The linear mixed-effects model estimating the random
                  effect showed some error messages/warnings. 3.Model misspecification.\n")
      }
     }else{resultid <- NULL} #Single individual will not have resultid table

    # Extract the excitation coeff for each excitation
    intdata[, totalexc := 0]
    intdata[, totalexcroll := 0]
    if(!noinput){
        if(verbose){print("Status: Calculating gains for one or several excitations")}
        for (i in 1:length(input)){ #For loop to go through all the inputs
          resultmean[, paste0(doremiexc[i],"_komega2") := summary$coefficients[paste0(doremiexc[i], "_rollmean"), "Estimate"]]
          resultmean[omega2>0, paste0(doremiexc[i],"_k") := summary$coefficients[paste0(doremiexc[i], "_rollmean"), "Estimate"]/resultmean$omega2]
          #If variation of the excitation coefficient across individuals needed:
          #And for each individual: the mean coeff (sumary$coeff) + the variation per Individual (in random)
          #Excitation coefficient in resultid
          if(nind > 1){  #Several individuals
            if(verbose){print("Status: several individuals")}

            resultid[, paste0(doremiexc[i],"_komega2") := (summary$coefficients[paste0(doremiexc[i], "_rollmean"), "Estimate"] +
                                                       random$id[.GRP,paste0(doremiexc[i], "_rollmean")]), by = id]
            resultid[omega2>0, paste0(doremiexc[i],"_k") := (summary$coefficients[paste0(doremiexc[i], "_rollmean"), "Estimate"] +
                                                                 random$id[.GRP,paste0(doremiexc[i], "_rollmean")])/omega2, by=id]

            intdata[, totalexc := totalexc +  resultid[[.GRP,paste0(doremiexc[i],"_k")]] * get(doremiexc[i]), by = id] #total excitation is the sum of all the excitations
            intdata[, totalexcroll := totalexcroll +  resultid[[.GRP,paste0(doremiexc[i],"_k")]] *  get(paste0(doremiexc[i],"_rollmean")), by = id] #total excitation is the sum of all the excitations

            #ponderated with their corresponding gains
          }else{ #Single individual
            if(verbose){print("Status: single individuals")}
            intdata[, totalexc := totalexc + resultmean[[paste0(doremiexc[i],"_k")]] * get(doremiexc[i])]
            intdata[, totalexcroll := totalexcroll + resultmean[[paste0(doremiexc[i],"_k")]] * get(paste0(doremiexc[i],"_rollmean"))]

          }
        }
    }else{ #there is no input
      if(verbose){print("Status: no input. Setting gain to 0.")}
      resultmean[, Komega2 := 0]
      if(nind>1){resultid[, Komega2 := 0]}
    }

    #The estimated signal is calculated by calling ode function in deSolve (through function "generate.1order"). As we will have a decomposition
    #of k for each excitation, the excitation considered is already the total excitation with the total gain (to avoid calculating both separately, this is why
    #k=1, total gain is already included in totalexc)
    #Assuming initial value is equilibrium value

    generate_solution<- function(time,time_derivate,totalexc,totalexcroll,signal_rollmean,signal_derivate1,xi,period,yeq){
      if (!(is.na(xi) | is.na(period) | is.na(yeq))){
        #Building time vector that had the first point of rollmean time and rollmean signal
        timecomp<-c(time[time<time_derivate[1]],time_derivate[1],time[time>time_derivate[1]])
        exccomp<-c(totalexc[time<time_derivate[1]],totalexcroll[1],totalexc[time>time_derivate[1]])

        #Position of the rollmean time in the original time vector
        i <- length(time[time<time_derivate[1]]) + 1
        droite <- generate.2order(time = timecomp[i:length(timecomp)], #includes original time and the first element of the rollmean time
                                  excitation = exccomp[i:length(exccomp)],
                                  y0 = signal_rollmean[1],
                                  v0 = signal_derivate1[1],
                                  xi = xi,
                                  period = period,
                                  k = 1,
                                  yeq = yeq)
        if(!dermethod == "calculate.fda"){ #with fda the initial condition will be on t=0 as we are building functions
          gauche <- generate.2order(time = timecomp[i:1],
                                  excitation = exccomp[i:1],
                                  y0 = signal_rollmean[1],
                                  v0 = signal_derivate1[1],
                                  xi = xi,
                                  period = period,
                                  k = 1,
                                  yeq = yeq)
        }else{gauche <-NULL}
        #putting the two solutions together
        sol<-c(gauche$y[i:2],droite$y)
       }else{
        if(verbose){print("Status: Solution wasn't generated as one of the parameters was NaN")}
        sol<-rep(NaN,length(time))
       }
      sol
    }
    if(nind > 1){
      if(verbose){print("Status: estimating signal for several individuals")}
      intdata[, signal_estimated := generate_solution(time,time_derivate,totalexc,totalexcroll,signal_rollmean,signal_derivate1,resultid[.GRP, xi],
                                          resultid[.GRP, period], resultid[.GRP, yeq]),by =id]

    }else{
      if(verbose){print("Status: estimating signal for single individual")}
      intdata[, signal_estimated := generate_solution(time,time_derivate,totalexc,totalexcroll,signal_rollmean,signal_derivate1,resultmean[, xi],
                                          resultmean[, period], resultmean[, yeq])]
    }

    }else{ # if the regression didn't work, a warning will be generated and tables will be set to NULL
    if (verbose){print("Status: Linear mixed-effect model produced errors.")}
    warning("Linear mixed-effect regression produced an error. Verify the regression object of the result.\n")
    resultid <- NULL
    resultmean <- NULL
    regression <- model
  }
  #Calculating R2
  intdata[, R2:= 1 - sum((signal - signal_estimated)^2,na.rm = T)/sum((signal - mean(signal,na.rm = T))^2,na.rm = T)]

  #Renaming columns in $data, $resultid, $resultmean objects to original names
  intdata[, id := NULL]
  if(!is.null(resultid)){resultid[, id := NULL]}

  #Replacing names in $data
  intdatanames <- names(intdata)
  intdatanamesnew <- names(intdata)

  rmeannames <- names(resultmean)
  rmeannamesnew <- names(resultmean)

  ridnames <- names(resultid)
  ridnamesnew <- names(resultid)

  for(idx in seq(doreminames)){
    intdatanamesnew <-gsub(paste0("^",doreminames[idx],"(?![0-9])"), originalnames[idx], intdatanamesnew, perl = T)
    #Replacing names in $resultmean
    if(!is.null(resultid)){
      ridnamesnew <-gsub(paste0("^",doreminames[idx],"(?![0-9])"), originalnames[idx], ridnamesnew, perl = T)
    }
    if(!is.null(resultmean)){
      rmeannamesnew <-gsub(paste0("^",doreminames[idx],"(?![0-9])"), originalnames[idx], rmeannamesnew, perl = T)
    }
  }
  setnames(intdata, intdatanames, intdatanamesnew)
  if(!is.null(resultid)){setnames(resultid, ridnames, ridnamesnew)}
  if(!is.null(resultmean)){setnames(resultmean, rmeannames, rmeannamesnew)}

  #Output the results of the function
  #Excitation string
  #If there is no excitation term
  if (noinput){str_exc <- 0
  }else{str_exc <- input} # If there is one OR SEVERAL excitation columns
  res = list(data = intdata, resultid = resultid, resultmean = resultmean, regression = regression, dermethod = dermethod, derparam = derparam, str_time = time, str_exc = str_exc, str_signal = signal, str_id = id)
  class(res)= "doremi" # Class definition
  return(res)
}

# Optimization loop for embedding dimension/spar --------------------------
#
optimum_param <- function(data,
                          id,
                          input= NULL,
                          time,
                          signal,
                          model = "analyze.1order",
                          dermethod = "calculate.gold", #string containing function name, beware!
                          pmin = 3,
                          pmax = 21,
                          pstep = 2,
                          verbose= F){
  analyze <- get(model)
  analysis <- rbindlist(lapply(seq(pmin,pmax,pstep),function(embedding){
    if(verbose){print(paste0("Analyzing for embedding=",embedding))}
    res <- analyze(data = data,
                      id = id,
                      input = input,
                      time = time,
                      signal = signal,
                      dermethod = dermethod,
                      derparam = embedding)
    res2 <- res$data
    res2[, xi_est:= res$resultmean$xi]
    res2[, period_est:= res$resultmean$period]
    res2[, D:= embedding]
    res2
  }))

  #Calculation of R2 max et Dopt
  R2opt <- max(analysis[!is.na(period_est),R2])
  Dopt <-analysis[R2==R2opt,D[1]]
  #Summary table for plotting
  summ<-analysis[,.(R2=R2[1],xi_est=xi_est[1],period_est=period_est[1]),by="D"]
  summ_opt<-summ[R2==R2opt]
  summ_opt[,method:=dermethod]
  names(summ_opt)[1]<-"Dopt"
  #Plotting estimated parameters and R2 versus embedding
  estvsembed<-ggplot(summ) + geom_point(aes(D,R2,color="R^2")) +
    geom_point(aes(D,xi_est,color="estimated damping factor, xi")) +
    geom_point(aes(D,period_est,color="estimated period, T"))+
    labs(x = "Embedding dimension, D",
         y = "",
         colour = "") +
    ggtitle("Evolution of estimated parameters and R^2 with embedding")+
    theme(legend.position = "top",
          plot.title = element_text(hjust = 0.5))
  return(list(analysis=analysis,summary=summ,summary_opt=summ_opt))
}



