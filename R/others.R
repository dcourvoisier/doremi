# Call to global variables ------------------------------------------------

globalVariables(c("id", "input", "time", "signal"))
globalVariables(c("id_tmp", "signal_derivate1", "signal_rollmean", "time_derivate"))
globalVariables(c("dampedsignalraw", "amplitudenorm", "dampedsignal", "excitation"))
globalVariables(c("exc_max", "exc_min", "intervals", "totalexc", "ymax", "ymin"))
globalVariables(c("tmpsignal", "expfit_A", "signal_estimated"))
globalVariables(c("value","variable"))
globalVariables(c("dampingtime", "eqvalue", "quantile", "sd","timedup"))

# Call to classes and methods ---------------------------------------------
# Classes documented directly in the analysis function (object "doremi")
# plot doremi generic method
#doremi <- function(x, ..., id) UseMethod("doremi")

#' S3 method to print DOREMI objects
#'
#' \code{print.doremi} prints the most important results of a DOREMI object
#' @param x DOREMI object
#' @param ... includes the additional arguments inherited from the generic print method
#' @return Returns the three coefficients of the differential equation estimated (fixed part, table $resultmean of the DOREMI object)
#' @examples
#' myresult <- remi(data = cardio,
#'                  id = "id",
#'                  input = "load",
#'                  time = "time",
#'                  signal = "hr",
#'                  embedding = 5)
#' myresult
#' @export
print.doremi = function (x, ...){
  print(x$resultmean)
}

#' S3 method to print DOREMI data objects
#'
#' \code{print.doremidata} prints the most important results of a  DOREMIDATA object
#' @param x DOREMIDATA object
#' @param ... includes the additional arguments inherited from the generic print method
#' @return Returns the table $data of the DOREMIDATA object
#' @examples
#' mydata <- generate.panel.remi(nind = 5,
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
  print("$data")
  print(x$data)
}

#' S3 method for DOREMI object summary
#'
#' \code{summary.doremi} provides a summary of the remi analysis
#' @param object DOREMI object (contains several lists)
#' @param ... includes the additional arguments inherited from the generic summary method
#' @return Returns a summary containing the five lists of the DOREMI object
#' @examples
#' myresult <- remi(data = cardio,
#'                  id = "id",
#'                  input = "load",
#'                  time = "time",
#'                  signal = "hr",
#'                  embedding = 5)
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
   cat("\n Estimated signal ($estimated):\n")
   print(object$estimated)
   cat("\n Embedding number ($embedding):\n")
   print(object$embedding)
}

#' S3 method to plot DOREMI objects
#'
#' \code{plot.doremi} generates a plot with the observed values of the signal, the excitation values and the fitted
#' signal over time for each individual.
#' @param x DOREMI object from \code{\link{remi}} analysis
#' @param ... includes the additional arguments inherited from the generic plot method
#' @param id Identifiers of the individuals to be represented in the plot.
#' By default, it will print the first six individuals.
#' @return Returns a plot with axis labels, legend and title. The axis labels and legend include the names of the variables set as input arguments.
#' The title includes the name of the DOREMI object result of the analysis. The function uses \code{\link[ggplot2]{ggplot}}
#' to generate the graphs and so it is possible to override the values of axis labels, legend and title through ggplot commands.
#' @examples
#' mydata <- generate.panel.remi(nind = 2,
#'                            dampingtime = 10,
#'                            amplitude = c(5,10),
#'                            nexc = 2,
#'                            duration = 20,
#'                            deltatf = 2,
#'                            tmax = 200,
#'                            minspacing = 0,
#'                            internoise = 0.2,
#'                            intranoise = 0.1)
#' myresult <- remi(data = mydata,
#'                  id = "id",
#'                  input = "excitation",
#'                  time = "time",
#'                  signal = "dampedsignal",
#'                  embedding = 5)
#' plot(myresult)
#' @export
#' @import ggplot2
#' @importFrom data.table melt
#' @importFrom data.table first
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

    if (first(x$str_exc != 0)){ #If there is an excitation term

      #Preparing a data frame in long format in order to use ggplot for all the excitations
      dataplot <- melt(xd, measure.vars = x$str_exc)

      #Excitations
      p <- p + geom_line(data = dataplot[get(x$str_id) %in% facets], aes(get(x$str_time), value, color = as.factor(variable)))
      #Estimated values
      if (is.null(xe[[x$str_id]])) {stop("id column not found in table \"estimated\". Check call to the analysis function, the id column
      should have been provided as input parameter.\n")}
      p <- p + geom_line(data = xe[get(x$str_id) %in% facets],
                          aes(get(x$str_time), ymin, colour = paste0(x$str_signal, " estimated, min value"))) +
               geom_line(data = xe[get(x$str_id) %in% facets],
                          aes(get(x$str_time), ymax, colour = paste0(x$str_signal, " estimated, max value")))

    }else{ #If there is no excitation term

      p <- p + geom_line(data = xe[get(x$str_id) %in% facets], aes(get(x$str_time), get(paste0(x$str_signal, "_estimated")), colour = paste0(x$str_signal, " estimated")))
    }

  }else{  # Single individual. Doremi object doesn't have id column

    p <- ggplot(xd)

    if (any(x$str_exc != 0)){ #If there is an excitation term
      #Preparing a data frame in long format in order to use ggplot for all the excitations
      dataplot <- melt(xd, measure.vars = x$str_exc)

      #Excitations
       p <- p + geom_line(data = dataplot, aes(get(x$str_time), value, color = as.factor(variable)))

      #Estimated values
       p <- p + geom_line(data = xe, aes(time, ymin, colour = paste0(x$str_signal, " estimated, min value"))) +
         geom_line(data = xe, aes(time, ymax, colour = paste0(x$str_signal, " estimated, max value")))

    }else{ #If there is no excitation term

      p <- p + geom_line(data = xe, aes(get(x$str_time), get(paste0(x$str_signal, "_estimated")), colour = paste0(x$str_signal, " estimated")))

    }

  }

  #Signal Can be 0 if using the function predict
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
#' \code{predict.doremi} predicts signal values with a DOREMI object when providing a new excitation vector(s).
#' @param object DOREMI object result of an analysis with the function remi
#' @param ... Additional arguments inherited from generic predict method.
#' @param newdata includes a data frame containing three columns or more:
#'
#' id (optional), indicating the individual identifier
#'
#'
#' time, containing the time values
#'
#' excitation, being one or several columns containing the different excitations
#' used to estimate a new signal. As in the other methods for the predict function, the columns of newdata
#' must have the same names as those of the original object.
#' @param verbose Is a boolean that displays status messages of the function when set to 1.
#' @return Returns a list containing the values of time, the values of the excitation and the predicted
#' values of the signal for the new excitation(s).
#' @examples
#' myresult <- remi(data = cardio[id == 1],
#'                  input = "load",
#'                  time = "time",
#'                  signal = "hr",
#'                  embedding = 5)
#' #Copying cardio into a new data frame and modifying the excitation column
#' new_exc <- cardio[id == 1, !"id"]
#'
#' et <- generate.excitation(amplitude = 100,
#'                                    nexc = 6,
#'                                    duration = 2,
#'                                    deltatf = 1,
#'                                    tmax = 49,
#'                                    minspacing = 2)
#' new_exc$load <- et$exc
#' new_exc$time <- et$t
#' predresult <- predict(myresult, newdata = new_exc)
#' plot(predresult)
#' @export
predict.doremi = function (object,
                           ...,
                           newdata,
                           verbose = FALSE){
  #Error management
  if (any(is.na(match(names(newdata), names(object$data))))) {
    stop("Cannot find the column names of data.frame 'newdata' in ", deparse(substitute(object))," .")
  }
  newdata <- setDT(newdata) #Convert data to data table

  if (!is.null(object$str_id)){ # Multiple individuals
      if (verbose){print("Predict status: panel data")}

      #Recovering id from object
      id <- object$str_id

      #Create column with time interval for each individual
      newdata[, deltat := 0.01 * min(diff(get(object$str_time))), by = id]

      #Calculation of the total excitation (sum of all the excitation columns times their regression coefficients)
      newdata[, totalexc := 0]
      for (i in seq(object$str_exc)){ #For loop to go through all the excitation columns
        newdata[, totalexc := totalexc + object$resultid[.GRP, get(paste0(object$str_exc[i], "_coeff"))]/
                    object$resultid[.GRP, dampingtime] * get(object$str_exc[i]), by = id]
      }

      # Generation of the third result table called \"estimated\"
      # Contains expanded time vector, minimum an maximum generated signals (for the two extreme scenarios of expanded excitation)
      #Generation of the expanded excitation vector according to desired deltat. Minimum and maximum values
      #Time vector and excitation vectors min and max
      #Expanded time vector: takes the minimum and the maximum of all the time intervals and creates the vector with deltat time intervals
      estimated <- newdata[, list(time = seq(floor(min(get(object$str_time), na.rm = T)), ceiling(max(get(object$str_time), na.rm = T)), deltat[1])), by = id]

      #Finding the values of the original time vector in the expanded time vector by using the function "findInterval"

      IDvec <- unique(newdata[[id]])

      for (idx in seq(IDvec)){ #It is necessary to do a loop because findInterval finds the index in which the value is found in the original time vector
        #And it will be necessary to shift these indexes according to in which id we are
        leftidx <- findInterval(newdata[get(object$str_id) == IDvec[idx], get(object$str_time)], estimated[get(object$str_id) == IDvec[idx], time])
        #tmpsignal contains the value of exctotal were the original times are found in the new time vector
        #and NA in the rest of the positions
        estimated[nrow(estimated[get(object$str_id) < IDvec[idx]]) + leftidx, tmpsignal := newdata[get(object$str_id) == IDvec[idx], totalexc]]

      }

      #The na.locf function (last observation called forward) will repeat the last non NA value.
      #We use it to repeat the values in the excitation function to the left and to the right
      estimated[, exc_min := na.locf(tmpsignal, na.rm = F, fromLast = F), by = id]
      estimated[, exc_max := na.locf(tmpsignal, na.rm = F, fromLast = T), by = id]

      #Calculating the convolution
      estimated[, ymin := generate.remi(object$resultid[.GRP, dampingtime],exc_min,time)$y +  object$resultid[.GRP, eqvalue], by = id]
      estimated[, ymax := generate.remi(object$resultid[.GRP, dampingtime],exc_max,time)$y +  object$resultid[.GRP, eqvalue], by = id]

  }else{  # Single individual
    if (verbose){print("Predict status: time series")}
    id <- NULL

    #Create column with time interval for each individual
    newdata[, deltat := 0.01 * min(diff(get(object$str_time)))]

    #Calculation of the total excitation (sum of all the excitation columns times their regression coefficients)
    newdata[, totalexc := 0]
    for (i in 1:length(object$str_exc)) #For loop to go through all the excitation columns
    {
      newdata[, totalexc := totalexc + object$resultmean[, get(paste0(object$str_exc[i], "_coeff"))] /
                                       object$resultmean[, dampingtime] * get(object$str_exc[i])]
    }

    # Expanded vectors according to detat chosen
    # Generation of the expanded excitation vector according to desired deltat. Minimum and maximum values
    # Time vector and excitation vectors min and max
    #Expanded time vector: takes the minimum and the maximum of all the time intervals and creates the vector with deltat time intervals
    estimated <- newdata[, list(time = seq(floor(min(get(object$str_time), na.rm = T)), ceiling(max(get(object$str_time), na.rm = T)), deltat[1]))]
    leftidx <- findInterval(newdata[, get(object$str_time)], estimated[, time])
    #tmpsignal contains the value of exctotal were the original times are found in the new time vector
    #and NA in the rest of the positions
    estimated[leftidx, tmpsignal := newdata[, totalexc]]

    #The na.locf function (last observation called forward) will repeat the last non NA value.
    #We use it to repeat the values in the excitation function to the left and to the right
    estimated[, exc_min := na.locf(tmpsignal, na.rm = F, fromLast = F)]
    estimated[, exc_max := na.locf(tmpsignal, na.rm = F, fromLast = T)]

    #Calculating the convolution
    estimated[, ymin := generate.remi(object$resultmean[, dampingtime],exc_min,time)$y + object$resultmean[, eqvalue]]
    estimated[, ymax := generate.remi(object$resultmean[, dampingtime],exc_max,time)$y + object$resultmean[, eqvalue]]
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


# Error management function -----------------------------------------------

#' Displays error messages for the analysis function according to the nature of the error
#' \code{errorcheck} displays error messages and/or warnings concerning the validity of input arguments provided to the analysis function
#' @param data data.frame or data.table containing the data to be analyzed. Same object that is passed as input argument to the analysis function.
#' @param col_var column variable. Contains a string that indicates the name of the column to analyze ("id","input",etc.)
#' @return Doesn't return a value. Either displays directly the error message/warning or changes data type in the data.frame/data.table provided
errorcheck = function (data, col_var){
  col_str <- deparse(substitute(col_var)) # Stores the name of the variable

  #It is a string but doesnt find a column with that name in the data table provided
  if(!(all(col_var %in% names(data))) & all(is.character(col_var))){
    stop("No column found in data with the parameter \"", col_str, "\" specified.\n")
  }
  #It is not a string
  if(!is.character(col_var)){
    if(col_str == "id"){contents <- "individual identifier"}
    if(col_str == "input"){contents <- "excitation"}
    if(col_str == "time"){contents <- "time"}
    if(col_str == "signal"){contents <- "signal"}
    stop(col_str," argument should be a string containing the name of the column in data that contains the ", contents,".\n")
  }
  #Verifies if data associated to the variable column is of the type string or factor and converts it to numeric
  #With the exception of id because for that case, the column is renamed and an extra column is created
  for (i in seq(col_var)){
    if(col_str!="id"){
      if (is.character(data[[col_var[i]]]) | is.factor(data[[col_var[i]]])){
        data[, col_var[i] := as.numeric(as.character(get(col_var[i])))]
        contents <- col_var[i]
        warning(contents," column was found to be of the factor/string type and was converted to numeric.\n")
      }
    }
  }
}
