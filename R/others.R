# Call to global variables ------------------------------------------------

globalVariables(c("inputcol", "timecol", "excitation_sum", "expfit_A"))
globalVariables(c("dampedsignalraw", "amplitudenorm", "dampedsignal", "exc_max", "exc_min", "excitation", "id", "intervals", "totalexc", "ymax", "ymin"))


# Call to classes and methods ---------------------------------------------
# Classes documented directly in the analysis function (object "doremi")
#' S3 method to plot DOREMI objects
#'
#' \code{plot.doremi} S3 method for the plot function so that it can represent DOREMI objects
#' @param doremi_object DOREMI object (contains several lists)
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
#' @importfrom graphics plot

plot.doremi = function (doremi_object,
                     id = 1:6){
  if (!is.null(doremi_object$data$id) && length(unique(doremi_object$data$id))>1){ # Multiple individuals
    facets <- as.character(id)
    p <- ggplot(doremi_object$data[doremi_object$data$id %in% facets]) +
         facet_wrap(~id, scales = "free")

    if (doremi_object$str_exc != 0){ #If there is an excitation term

      #Preparing a data frame in long format in order to use ggplot for all the excitations
      dataplot <- melt(doremi_object$data, measure.vars = doremi_object$str_exc)

      #Excitations
      p <- p + geom_line(data = dataplot[dataplot$id %in% facets], aes(get(doremi_object$str_time), value, color = as.factor(variable)))

      #Estimated values
      p <- p + geom_line(data = doremi_object$estimated[doremi_object$estimated$id %in% facets],
                           aes(timecol, ymin, colour = paste0(doremi_object$str_signal, " estimated, min value"))) +
                 geom_line(data = doremi_object$estimated[doremi_object$estimated$id %in% facets],
                           aes(timecol, ymax, colour = paste0(doremi_object$str_signal, " estimated, max value")))

    }else{ #If there is no excitation term

      p <- p + geom_line(data = doremi_object$estimated[doremi_object$estimated$id %in% facets], aes(get(doremi_object$str_time), get(paste0(doremi_object$str_signal, "_estimated")), colour = paste0(doremi_object$str_signal, " estimated")))
    }

  }else{  # Single individual. Doremi object doesn't have id column

    p <- ggplot(doremi_object$data)

    if (doremi_object$str_exc != 0){ #If there is an excitation term
      #Preparing a data frame in long format in order to use ggplot for all the excitations
      dataplot <- melt(doremi_object$data, measure.vars = doremi_object$str_exc)

      #Excitations
       p <- p + geom_line(data = dataplot, aes(get(doremi_object$str_time), value, color = as.factor(variable)))

      #Estimated values
       p <- p + geom_line(data = doremi_object$estimated, aes(timecol, ymin, colour = paste0(doremi_object$str_signal, " estimated, min value"))) +
         geom_line(data = doremi_object$estimated, aes(timecol, ymax, colour = paste0(doremi_object$str_signal, " estimated, max value")))

    }else{ #If there is no excitation term

      p <- p + geom_line(data = doremi_object$estimated, aes(get(doremi_object$str_time), get(paste0(doremi_object$str_signal, "_estimated")), colour = paste0(doremi_object$str_signal, " estimated")))

    }

  }

  #Signal. Can be 0 if using the function predict
  if (!is.null(doremi_object$str_signal)){
    p <- p + geom_point(aes(get(doremi_object$str_time), get(doremi_object$str_signal), colour = paste0("Signal: ", doremi_object$str_signal)))
  }

  p <- p + labs(x = doremi_object$str_time,
       y = doremi_object$str_signal,
       colour = "") +
  ggtitle(paste0("DOREMI plot: ",deparse(substitute(doremi_object))))+
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5))
  return(p)
}

#' S3 method to predict signal values in a DOREMI object when entering a new excitation
#'
#' \code{predict.doremi} S3 method for the predict function so that it can
#' predict signal values in a DOREMI object when entering a new excitation
#' @param doremi_object DOREMI object result of an analysis with the function doremi_analyse_order1
#' @param newdata data table containing three columns or more: id (optional), indicating the individual identifier,
#' time, containing the time values and excitation, being one or several columns containing the different excitations
#' to which the user desires to obtain an estimated signal. As in the other methods for the predict function, the columns of newdata
#' must have the same names as those of the original object, otherwise an error message is displayed
#' @return Returns an list containing the values of time, the values of the excitation and the predicted
#' values of the signal for that excitation
#' @examples
#' @export
#' @importfrom stats predict
#' @importfrom data.table .GRP
#' @importfrom data.table :=

predict.doremi = function (doremi_object,
                           newdata){
  deltat <- 0.1
  #Error management
  if (any(is.na(match(names(newdata), names(doremi_object$data))))) {
    stop("cannot evaluate groups for desired levels on 'newdata'")
  }

  newdata <- setDT(newdata)
  if (!is.null(newdata$id)){ # Multiple individuals

      #Calculation of the total excitation (sum of all the excitation columns times their regression coefficients)
      newdata[, totalexc := 0]
      for (i in 1:length(doremi_object$str_exc)){ #For loop to go through all the excitation columns
        newdata[, totalexc := totalexc + doremi_object$resultid[.GRP, get(paste0(doremi_object$str_exc[i], "_exccoeff"))]/
                    doremi_object$resultid[.GRP, get(paste0(doremi_object$str_signal, "_dampingtime"))] *
                    get(doremi_object$str_exc[i]), by = id]
      }
      # Generation of the third result table called \"estimated\"
      # Contains expanded time vector, minimum an maximum generated signals (for the two extreme scenarios of expanded excitation)
      #Generation of the expanded excitation vector according to desired deltat. Minimum and maximum values
      #Time vector and excitation vectors min and max
      #Expanded time vector: takes the minimum and the maximum of all the time intervals and creates the vector with deltat time intervals
      estimated <- newdata[, list(timecol = seq(floor(min(get(doremi_object$str_time), na.rm = T)), ceiling(max(get(doremi_object$str_time), na.rm = T)), deltat)), by = id]

      #Finding the values of the original time vector in the expanded time vector by using the function "findInterval"
      IDvec <- unique(newdata$id)

      for (idx in seq(IDvec)){ #It is necessary to do a loop because findInterval finds the index in which the value is found in the original time vector
        #And it will be necessary to shift these indexes according to in which id we are
        leftidx <- findInterval(newdata[id == IDvec[idx], get(doremi_object$str_time)], estimated[id == IDvec[idx], timecol])
        #tmpsignal contains the value of exctotal were the original times are found in the new time vector
        #and NA in the rest of the positions
        estimated[nrow(estimated[id < IDvec[idx]]) + leftidx, tmpsignal := newdata[id == IDvec[idx], totalexc]]

      }

      #The na.locf function (last observation called forward) will repeat the last non NA value.
      #We use it to repeat the values in the excitation function to the left and to the right
      estimated[, exc_min := na.locf(tmpsignal, na.rm = F, fromLast = F), by = id]
      estimated[, exc_max := na.locf(tmpsignal, na.rm = F, fromLast = T), by = id]

      #Calculating the convolution
      estimated[, ymin := doremi_generate_order1(doremi_object$resultid[.GRP, get(paste0(doremi_object$str_signal,"_dampingtime"))],exc_min,timecol)$y+
                  doremi_object$resultid[.GRP, get(paste0(doremi_object$str_signal,"_eqvalue"))], by = id]
      estimated[, ymax := doremi_generate_order1(doremi_object$resultid[.GRP, get(paste0(doremi_object$str_signal,"_dampingtime"))],exc_max,timecol)$y+
                  doremi_object$resultid[.GRP, get(paste0(doremi_object$str_signal,"_eqvalue"))], by = id]
  }else{  # Single individual

    #Calculation of the total excitation (sum of all the excitation columns times their regression coefficients)
    newdata[, totalexc := 0]
    for (i in 1:length(doremi_object$str_exc)) #For loop to go through all the excitation columns
    {
      newdata[, totalexc := totalexc + doremi_object$resultmean[, get(paste0(doremi_object$str_exc[i], "_exccoeff"))]/
                                         doremi_object$resultmean[, get(paste0(doremi_object$str_signal, "_dampingtime"))] *
                                         get(doremi_object$str_exc[i])]
    }

    # Expanded vectors according to detat chosen
    # Generation of the expanded excitation vector according to desired deltat. Minimum and maximum values
    # Time vector and excitation vectors min and max
    #Expanded time vector: takes the minimum and the maximum of all the time intervals and creates the vector with deltat time intervals
    estimated <- doremi_object$data[, list(timecol = seq(floor(min(get(doremi_object$str_time), na.rm = T)), ceiling(max(get(doremi_object$str_time), na.rm = T)), deltat))]
    leftidx <- findInterval(newdata[, get(doremi_object$str_time)], estimated[, timecol])
    #tmpsignal contains the value of exctotal were the original times are found in the new time vector
    #and NA in the rest of the positions
    estimated[leftidx, tmpsignal := newdata[, totalexc]]

    #The na.locf function (last observation called forward) will repeat the last non NA value.
    #We use it to repeat the values in the excitation function to the left and to the right
    estimated[, exc_min := na.locf(tmpsignal, na.rm = F, fromLast = F)]
    estimated[, exc_max := na.locf(tmpsignal, na.rm = F, fromLast = T)]

    #Calculating the convolution
    estimated[, ymin := doremi_generate_order1(doremi_object$resultmean[, get(paste0(doremi_object$str_signal, "_dampingtime"))],exc_min,timecol)$y+
                doremi_object$resultmean[, get(paste0(doremi_object$str_signal,"_eqvalue"))]]
    estimated[, ymax := doremi_generate_order1(doremi_object$resultmean[, get(paste0(doremi_object$str_signal, "_dampingtime"))],exc_max,timecol)$y+
                doremi_object$resultmean[, get(paste0(doremi_object$str_signal,"_eqvalue"))]]

  }
  #Creating a DOREMI object to return so that it can be easily plotted afterwards
  #This object will be a copy of the original doremi_object set as input,
  #however the following attributes will be updated:
  #data - will now be reduced to the inputs given to the function: newdata and the added column totalexc
  #estimated - will now contain the estimated signal min and max values for the excitation introduced as input
  #str_time will be updated with name of time column in table newdata
  #str_exc will be updated with name of excitation columns in table newdata

  res <- doremi_object
  res$data <- newdata
  res$estimated <- estimated
  res$str_time <- doremi_object$str_time
  res$str_exc <- doremi_object$str_exc
  res$str_signal <- NULL


  return(res)
}
