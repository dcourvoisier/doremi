# Call to global variables ------------------------------------------------

globalVariables(c("inputcol", "timecol", "excitation_sum", "expfit_A"))
globalVariables(c("dampedsignalraw", "amplitudenorm", "dampedsignal", "exc_max", "exc_min", "excitation", "id", "intervals", "totalexc", "ymax", "ymin"))
globalVariables(c("idchar","tmpsignal"))
globalVariables(c("value","variable"))

# Call to classes and methods ---------------------------------------------
# Classes documented directly in the analysis function (object "doremi")
# plot doremi generic method
#doremi <- function(x, ..., id) UseMethod("doremi")

#' S3 method to print DOREMI objects
#'
#' \code{print.doremi} S3 method for the print function so that it can represent DOREMI objects
#' @param x DOREMI object (contains several lists)
#' @param ... includes the additional arguments inherited from the generic print method
#' @return Returns the main parameters of the DOREMI object: the mean values of the three coefficients of the differential equation
#' @examples
#' mydata <- simulation_generate_order1(nind = 5,
#'                            dampingtime = 10,
#'                            amplitude = c(5,10),
#'                            nexc = 2,
#'                            duration = 20,
#'                            deltatf = 2,
#'                            tmax = 200,
#'                            minspacing = 0,
#'                            internoise = 0.2,
#'                            intranoise = 0.1)
#' myresult <- doremi_analyse_order1(data = mydata,
#'                            id = "id",
#'                            input = "excitation",
#'                            time = "timecol",
#'                            signal = "dampedsignal",
#'                            embedding = 5)
#' myresult
#' @export
print.doremi = function (x, ...){
  print(x$resultmean)
}

#' S3 method to print DOREMI data objects
#'
#' \code{print.doremidata} S3 method for the print function so that it can show DOREMI data objects
#' @param x DOREMI object (contains several lists)
#' @param ... includes the additional arguments inherited from the generic print method
#' @return Returns the main parameters of the DOREMI data object
#' @examples
#' mydata <- simulation_generate_order1(nind = 5,
#'                            dampingtime = 10,
#'                            amplitude = c(5,10),
#'                            nexc = 2,
#'                            duration = 20,
#'                            deltatf = 2,
#'                            tmax = 200,
#'                            minspacing = 0,
#'                            internoise = 0.2,
#'                            intranoise = 0.1)
#' mydata
#' @export
print.doremidata = function (x, ...){
  print(x$data)
}

#' S3 method for DOREMI object summary
#'
#' \code{summary.doremi} S3 method for the summary function so that it can represent DOREMI objects
#' @param object DOREMI object (contains several lists)
#' @param ... includes the additional arguments inherited from the generic summary method
#' @return Returns a summary with all the lists of the DOREMI object
#' @examples
#' mydata <- simulation_generate_order1(nind = 5,
#'                            dampingtime = 10,
#'                            amplitude = c(5,10),
#'                            nexc = 2,
#'                            duration = 20,
#'                            deltatf = 2,
#'                            tmax = 200,
#'                            minspacing = 0,
#'                            internoise = 0.2,
#'                            intranoise = 0.1)
#' myresult <- doremi_analyse_order1(data = mydata,
#'                            id = "id",
#'                            input = "excitation",
#'                            time = "timecol",
#'                            signal = "dampedsignal",
#'                            embedding = 5)
#' summary(myresult)
#' @export
summary.doremi = function (object, ...){
   cat("Derivative and mean calculation ($data):\n")
   print(object$data)
   cat("\n Mixed-effects regression results ($regression):\n")
   print(object$regression)
   cat("\n Mean coefficients of the differential equation ($resultmean):\n")
   print(object$resultmean)
   cat("\n Coefficients per individual ($resultid):\n")
   print(object$resultid)
}

#' S3 method to plot DOREMI objects
#'
#' \code{plot.doremi} S3 method for the plot function so that it can represent DOREMI objects
#' @param x DOREMI object (contains several lists)
#' @param ... includes the additional arguments inherited from the generic plot method
#' @param id Includes the identifiers of the individuals to be represented in the plot.
#' id can be a scalar (the function will plot an individual) or a vector (it will plot the individuals with
#' the numeric id contained in the vector)
#' By default, it will print the six first individuals.
#' @return Returns a plot with legend and title that includes the name of the DOREMI object
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
#'                            signal = "dampedsignal",
#'                            embedding = 5)
#' plot(myresult)
#' plot(myresult,1)
#' plot(myresult,2:5)
#' @export
#' @import ggplot2
#' @importFrom data.table melt
plot.doremi = function (x, ...,
                     id = NULL){
  xd <- x$data # data frame that contains the data from the doremi result object
  xe <- x$estimated # data frame containing the signal estimated values from the doremi result object

  if (!is.null(x$str_id) && length(unique(xd[[x$str_id]])) > 1){ # Multiple individuals
    if (is.null(id)){#if user doesn't specify id, plot prints the first six individuals
      facets <- head(unique(xd[[x$str_id]]), 6)
    }else{
      facets <- as.character(id)
    }
    p <- ggplot(xd[get(x$str_id) %in% facets]) +
         facet_wrap(~get(x$str_id), scales = "free")

    if (x$str_exc != 0){ #If there is an excitation term

      #Preparing a data frame in long format in order to use ggplot for all the excitations
      dataplot <- melt(xd, measure.vars = x$str_exc)

      #Excitations
      p <- p + geom_line(data = dataplot[get(x$str_id) %in% facets], aes(get(x$str_time), value, color = as.factor(variable)))

      #Estimated values
      if (is.null(xe[[x$str_id]])) {stop("id column not found in table \"estimated\". Check call to the analysis function, the id column
      should have been provided as input parameter.\n")}
      p <- p + geom_line(data = xe[get(x$str_id) %in% facets],
                          aes(timecol, ymin, colour = paste0(x$str_signal, " estimated, min value"))) +
               geom_line(data = xe[get(x$str_id) %in% facets],
                          aes(timecol, ymax, colour = paste0(x$str_signal, " estimated, max value")))

    }else{ #If there is no excitation term

      p <- p + geom_line(data = xe[get(x$str_id) %in% facets], aes(get(x$str_time), get(paste0(x$str_signal, "_estimated")), colour = paste0(x$str_signal, " estimated")))
    }

  }else{  # Single individual. Doremi object doesn't have id column

    p <- ggplot(xd)

    if (x$str_exc != 0){ #If there is an excitation term
      #Preparing a data frame in long format in order to use ggplot for all the excitations
      dataplot <- melt(xd, measure.vars = x$str_exc)

      #Excitations
       p <- p + geom_line(data = dataplot, aes(get(x$str_time), value, color = as.factor(variable)))

      #Estimated values
       p <- p + geom_line(data = xe, aes(timecol, ymin, colour = paste0(x$str_signal, " estimated, min value"))) +
         geom_line(data = xe, aes(timecol, ymax, colour = paste0(x$str_signal, " estimated, max value")))

    }else{ #If there is no excitation term

      p <- p + geom_line(data = xe, aes(get(x$str_time), get(paste0(x$str_signal, "_estimated")), colour = paste0(x$str_signal, " estimated")))

    }

  }

  #Signal. Can be 0 if using the function predict
  if (!is.null(x$str_signal)){
    p <- p + geom_point(aes(get(x$str_time), get(x$str_signal), colour = x$str_signal))
  }

  p <- p + labs(x = x$str_time,
       y = x$str_signal,
       colour = "") +
  ggtitle(paste0("DOREMI plot: ",deparse(substitute(x))))+
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5))
  return(p)
}

#' S3 method to predict signal values in a DOREMI object when entering a new excitation
#'
#' \code{predict.doremi} S3 method for the predict function so that it can
#' predict signal values in a DOREMI object when entering a new excitation
#' @param object DOREMI object result of an analysis with the function doremi_analyse_order1
#' @param ... Additionnal arguments inherited from generic predict method.
#' @param newdata includes a data table containing three columns or more: id (optional), indicating the individual identifier,
#' time, containing the time values and excitation, being one or several columns containing the different excitations
#' to which the user desires to obtain an estimated signal. As in the other methods for the predict function, the columns of newdata
#' must have the same names as those of the original object, otherwise an error message is displayed
#' @return Returns an list containing the values of time, the values of the excitation and the predicted
#' values of the signal for that excitation
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
#'                            signal = "dampedsignal",
#'                            embedding = 5)
#' exc_table <- myresult$data
#' predresult <- predict(myresult, newdata = exc_table)
#' plot(predresult)
#' @export
predict.doremi = function (object, ..., newdata){
  deltat <- 0.1
  #Error management
  if (any(is.na(match(names(newdata), names(object$data))))) {
    stop("Cannot evaluate groups for desired levels on 'newdata'")
  }

  newdata <- setDT(newdata) #COnvert data to data table
  id <- object$str_id #Recover id name from the object str_id attribute

  if (!is.null(id)){ # Multiple individuals

      #Calculation of the total excitation (sum of all the excitation columns times their regression coefficients)
      newdata[, totalexc := 0]
      for (i in 1:length(object$str_exc)){ #For loop to go through all the excitation columns
        newdata[, totalexc := totalexc + object$resultid[.GRP, get(paste0(object$str_exc[i], "_exccoeff"))]/
                    object$resultid[.GRP, get(paste0(object$str_signal, "_dampingtime"))] *
                    get(object$str_exc[i]), by = get(object$str_id)]
      }
      # Generation of the third result table called \"estimated\"
      # Contains expanded time vector, minimum an maximum generated signals (for the two extreme scenarios of expanded excitation)
      #Generation of the expanded excitation vector according to desired deltat. Minimum and maximum values
      #Time vector and excitation vectors min and max
      #Expanded time vector: takes the minimum and the maximum of all the time intervals and creates the vector with deltat time intervals
      estimated <- newdata[, list(timecol = seq(floor(min(get(object$str_time), na.rm = T)), ceiling(max(get(object$str_time), na.rm = T)), deltat)), by = id]

      #Finding the values of the original time vector in the expanded time vector by using the function "findInterval"

      IDvec <- unique(newdata[[id]])

      for (idx in seq(IDvec)){ #It is necessary to do a loop because findInterval finds the index in which the value is found in the original time vector
        #And it will be necessary to shift these indexes according to in which id we are
        leftidx <- findInterval(newdata[get(object$str_id) == IDvec[idx], get(object$str_time)], estimated[get(object$str_id) == IDvec[idx], timecol])
        #tmpsignal contains the value of exctotal were the original times are found in the new time vector
        #and NA in the rest of the positions
        estimated[nrow(estimated[get(object$str_id) < IDvec[idx]]) + leftidx, tmpsignal := newdata[get(object$str_id) == IDvec[idx], totalexc]]

      }
      print(estimated)
      #The na.locf function (last observation called forward) will repeat the last non NA value.
      #We use it to repeat the values in the excitation function to the left and to the right
      estimated[, exc_min := na.locf(tmpsignal, na.rm = F, fromLast = F), by = id]
      estimated[, exc_max := na.locf(tmpsignal, na.rm = F, fromLast = T), by = id]

      #Calculating the convolution
      estimated[, ymin := doremi_generate_order1(object$resultid[.GRP, get(paste0(object$str_signal,"_dampingtime"))],exc_min,timecol)$y+
                  object$resultid[.GRP, get(paste0(object$str_signal,"_eqvalue"))], by = id]
      estimated[, ymax := doremi_generate_order1(object$resultid[.GRP, get(paste0(object$str_signal,"_dampingtime"))],exc_max,timecol)$y+
                  object$resultid[.GRP, get(paste0(object$str_signal,"_eqvalue"))], by = id]
  }else{  # Single individual

    #Calculation of the total excitation (sum of all the excitation columns times their regression coefficients)
    newdata[, totalexc := 0]
    for (i in 1:length(object$str_exc)) #For loop to go through all the excitation columns
    {
      newdata[, totalexc := totalexc + object$resultmean[, get(paste0(object$str_exc[i], "_exccoeff"))]/
                                         object$resultmean[, get(paste0(object$str_signal, "_dampingtime"))] *
                                         get(object$str_exc[i])]
    }

    # Expanded vectors according to detat chosen
    # Generation of the expanded excitation vector according to desired deltat. Minimum and maximum values
    # Time vector and excitation vectors min and max
    #Expanded time vector: takes the minimum and the maximum of all the time intervals and creates the vector with deltat time intervals
    estimated <- object$data[, list(timecol = seq(floor(min(get(object$str_time), na.rm = T)), ceiling(max(get(object$str_time), na.rm = T)), deltat))]
    leftidx <- findInterval(newdata[, get(object$str_time)], estimated[, timecol])
    #tmpsignal contains the value of exctotal were the original times are found in the new time vector
    #and NA in the rest of the positions
    estimated[leftidx, tmpsignal := newdata[, totalexc]]

    #The na.locf function (last observation called forward) will repeat the last non NA value.
    #We use it to repeat the values in the excitation function to the left and to the right
    estimated[, exc_min := na.locf(tmpsignal, na.rm = F, fromLast = F)]
    estimated[, exc_max := na.locf(tmpsignal, na.rm = F, fromLast = T)]

    #Calculating the convolution
    estimated[, ymin := doremi_generate_order1(object$resultmean[, get(paste0(object$str_signal, "_dampingtime"))],exc_min,timecol)$y+
                object$resultmean[, get(paste0(object$str_signal,"_eqvalue"))]]
    estimated[, ymax := doremi_generate_order1(object$resultmean[, get(paste0(object$str_signal, "_dampingtime"))],exc_max,timecol)$y+
                object$resultmean[, get(paste0(object$str_signal,"_eqvalue"))]]

  }
  #Creating a DOREMI object to return so that it can be easily plotted afterwards
  #This object will be a copy of the original object set as input,
  #however the following attributes will be updated:
  #data - will now be reduced to the inputs given to the function: newdata and the added column totalexc
  #estimated - will now contain the estimated signal min and max values for the excitation introduced as input
  #str_time will be updated with name of time column in table newdata
  #str_exc will be updated with name of excitation columns in table newdata

  res <- object
  res$data <- newdata
  res$estimated <- estimated
  res$str_time <- object$str_time
  res$str_exc <- object$str_exc
  res$str_signal <- NULL
  res$str_id <- id


  return(res)
}

