# Call to global variables ------------------------------------------------

globalVariables(c("id", "input", "time", "signal","signalraw","excitation"))
globalVariables(c("id_tmp", "signal_derivate1", "signal_rollmean", "time_derivate"))
globalVariables(c("totalexc", "totalexcroll"))
globalVariables(c("tmpsignal", "expfit_A", "signal_estimated","tau","timedup","value","variable"))
globalVariables(c("Komega2","R2","deltatf","esp2omega","method","omega2","period", "signal_derivate2","wn","xi","yeq","yeqomega2"))

# Call to classes and methods ---------------------------------------------

#' S3 method to print DOREMI objects
#'
#' \code{print.doremi} prints the most important results of a DOREMI object
#' @param x DOREMI object
#' @param ... includes the additional arguments inherited from the generic print method
#' @return Returns the coefficients of the differential equation estimated (fixed coefficients, table $resultmean of the DOREMI object)
#' @examples
#' myresult <- analyze.1order(data = cardio,
#'                  id = "id",
#'                  input = "load",
#'                  time = "time",
#'                  signal = "hr")
#' myresult
#' @export
print.doremi = function (x, ...){
  print(x$resultmean)
}

#' S3 method to print DOREMIDATA objects
#'
#' \code{print.doremidata} prints the a DOREMIDATA object
#' @param x DOREMIDATA object
#' @param ... includes the additional arguments inherited from the generic print method
#' @return Returns the DOREMIDATA object (datatable))
#' @examples
#' time <- 0:100
#' data <- generate.panel.2order(time = time,
#'                               y0 = 10,
#'                               v0 = 0,
#'                               xi = 0.1,
#'                               period = 30,
#'                               k = 1,
#'                               yeq = 2,
#'                               nind = 6,
#'                               internoise = 0.3,
#'                               intranoise = 5)
#' data
#' @export
print.doremidata = function (x, ...){
  class(x)<-"data.table"
  print(x)
}
#' S3 method for DOREMI object summary
#'
#' \code{summary.doremi} provides a summary containing the five lists of the DOREMI object
#' @param object, DOREMI object (contains several lists)
#' @param ... includes the additional arguments inherited from the generic summary method
#' @return Returns a summary containing the five lists of the DOREMI object
#' @examples
#' myresult <- analyze.1order(data = cardio,
#'                  id = "id",
#'                  input = "load",
#'                  time = "time",
#'                  signal = "hr")
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
   cat("\n Derivative estimation method used ($dermethod):\n")
   print(object$dermethod)
   cat("\n Embedding number/smoothing parameter ($derparam):\n")
   print(object$derparam)
}

#' S3 method to plot DOREMI objects
#'
#' \code{plot.doremi} generates a plot with the observed values of the signal, the excitation values and the fitted
#' signal over time for each individual.
#' @param x, DOREMI object resulting from \code{\link{analyze.1order}} or \code{\link{analyze.2order}} analysis
#' @param ... includes the additional arguments inherited from the generic plot method
#' @param id Identifiers of the individuals to be represented in the plot.
#' By default, it will print the first six individuals.
#' @return Returns a plot with axis labels, legend and title. The axis labels and legend include the names of the variables set as input arguments.
#' The title includes the name of the DOREMI object result of the analysis. The function uses \code{\link[ggplot2]{ggplot}}
#' to generate the graphs and so it is possible to override the values of axis labels, legend and title through ggplot commands.
#' @examples
#' mydata <- generate.panel.1order(time= 0:100,
#'                                 excitation = sin(0:100),
#'                                 y0 = 0,
#'                                 t0 = 0,
#'                                 tau = 2,
#'                                 k = 1,
#'                                 yeq = 0,
#'                                 nind = 2,
#'                                 internoise = 0.1,
#'                                 intranoise = 8)
#' myresult <- analyze.1order(data = mydata,
#'                            id = "id",
#'                            input = "excitation",
#'                            time = "time",
#'                            signal = "signal")
#' plot(myresult)
#' @export
#' @import ggplot2
#' @importFrom data.table melt
#' @importFrom data.table first
plot.doremi = function (x, ...,
                     id = NULL){
  xd <- x$data # data frame that contains the data from the doremi result object

  if (!is.null(x$str_id) && length(unique(xd[[x$str_id]])) > 1){ # Multiple individuals
    if (is.null(id)){#if user doesn't specify id, plot prints the first six individuals
      facets <- head(unique(xd[[x$str_id]]), 6)
    }else{
      facets <- as.character(id)
    }
  }else{#Single individual
      facets <- 1}

    p <- ggplot(xd[get(x$str_id) %in% facets]) +
         facet_wrap(~get(x$str_id), scales = "free")

    if (first(x$str_exc != 0)){ #If there is an excitation term
      #Preparing a data frame in long format in order to use ggplot for all the excitations
      dataplot <- melt(xd, measure.vars = x$str_exc)
      #Excitations
      p <- p + geom_line(data = dataplot[get(x$str_id) %in% facets], aes(get(x$str_time), value, color = as.factor(variable)))
    }

    #Estimated values
    p <- p + geom_line(aes(get(x$str_time), get(paste0(x$str_signal,"_estimated")), colour = paste0(x$str_signal, " estimated")))

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
#' S3 method to plot DOREMIDATA objects
#'
#' \code{plot.doremidata} generates a plot of the simulated signals resulting from the \code{\link{generate.panel.1order}} and \code{\link{generate.panel.2order}} functions
#' @param x DOREMIDATA object resulting from the aforementioned functions
#' @param ... includes the additional arguments inherited from the generic plot method
#' @return Returns a plot with axis labels, legend and title. The title includes the name of the DOREMIDATA object result of the analysis.
#' The function uses \code{\link[ggplot2]{ggplot}}
#' to generate the graphs and thus it is possible to override the values of axis labels, legend and title through ggplot commands.
#' @examples
#' mydata <- generate.panel.1order(time=0:100,
#'                                 excitation = c(rep(0,50),rep(1,51)),
#'                                 nind = 6,
#'                                 internoise = 0.2,
#'                                 intranoise = 100)
#' plot(mydata)
#' @export
#' @import ggplot2
plot.doremidata = function (x, ...){

  p <- ggplot(x) + geom_point(aes(time, signal, color = "signal")) +
    geom_line(aes(time, signalraw, color = "signal, no noise")) +
    facet_wrap(~id, scales = "free") +
    labs(x = "time",
                y = "signal",
                colour = "") +
    ggtitle(paste0("DOREMI simulated data: ",deparse(substitute(x))))+
    theme(legend.position = "top",
          plot.title = element_text(hjust = 0.5))
  if(!is.null(x$excitation)){
    p <- p + geom_line(aes(time, excitation, color = "excitation"))
  }
  return(p)

}
#' S3 method to plot DOREMIPARAM objects
#'
#' \code{plot.doremiparam} generates a plot of the parameters resulting from the \code{\link{optimum_param}} function
#' @param x DOREMIPARAM object resulting from the aforementioned function
#' @param ... includes the additional arguments inherited from the generic plot method
#' @return Returns a plot showing the evolution of the first/second order differential equation coefficients and R2 with the values taken by the embedding number/smoothing parameter
#' (see details of \code{\link{optimum_param}} function).
#' The function uses \code{\link[ggplot2]{ggplot}}
#' to generate the graphs and thus it is possible to override the values of axis labels, legend and title through ggplot commands.
#' @examples
#' mydata <- generate.panel.1order(time = 0:130,
#'                           excitation = c(rep(0,30),rep(1,50),rep(0,51)),
#'                           nind = 5,
#'                           internoise = 0.2,
#'                           intranoise = 100)
#' myres<- optimum_param (data = mydata,
#'                          id = "id",
#'                          input ="excitation",
#'                          time = "time",
#'                          signal = "signal",
#'                          model = "1order",
#'                          dermethod = "gold",
#'                          pmin = 3,
#'                          pmax = 11,
#'                          pstep = 2)
#' plot(myres)
#' @export
#' @import ggplot2
#' @importFrom data.table melt
plot.doremiparam = function (x, ...){
  analysis <- x$analysis
  dermethod <-x$summary_opt$method
  #Temporary long data table containing parameter name
  toplot<-melt(analysis[,-c("id")],id.vars="D")
  #Plotting estimated parameters and R2 versus embedding
  estvsembed<-ggplot(toplot) +
    geom_point(aes(D,value,color=variable)) +
    labs(x = "Embedding dimension, D",
         y = "",
         colour = "") +
    ggtitle(paste0("Evolution of R2 and the estimated parameters\nwith the embedding dimension/smoothing parameter: ",dermethod))+
    theme_bw()+
    theme(legend.position = "top",
          plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~variable,scales = "free")
  return(estvsembed)
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
#' myresult <- analyze.1order(data = cardio[id==1],
#'                  id="id",
#'                  input = "load",
#'                  time = "time",
#'                  signal = "hr")
#' #Copying cardio into a new data frame and modifying the excitation column
#' new_exc <- cardio[id==1]
#' et <- generate.excitation(amplitude = 100,
#'                           nexc = 6,
#'                           duration = 2,
#'                           deltatf = 1,
#'                           tmax = 49,
#'                           minspacing = 2)
#' new_exc$load <- et$exc
#' new_exc$time <- et$t
#' predresult <- predict(myresult, newdata = new_exc)
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
  #Initialization of the total excitation (sum of all the excitation columns times their regression coefficients)
  newdata[, totalexc := 0]

  if (!is.null(object$str_id) & length(unique(object$data[[object$str_id]])) > 1){ # Multiple individuals
      if (verbose){print("Predict status: panel data")}
      #Recovering id from object
      id <- object$str_id
      for (i in seq(object$str_exc)){ #For loop to go through all the excitation columns
        newdata[, totalexc := totalexc + object$resultid[.GRP, get(paste0(object$str_exc[i], "_k"))] * get(object$str_exc[i]), by = id]
      }
      #Calculating the predicted signal
      newdata[, paste0(object$str_signal, "_estimated") := generate.1order(time = get(object$str_time),
                                                                           excitation = totalexc,
                                                                           y0 = object$resultid[.GRP, yeq],
                                                                           tau = object$resultid[.GRP, tau],
                                                                           k = 1,
                                                                           yeq = object$resultid[.GRP, yeq])$y, by = id]
  }else{  # Single individual
      id <- "id"
      if (verbose){print("Predict status: time series")}
      for (i in 1:length(object$str_exc)){ #For loop to go through all the excitation columns
        newdata[, totalexc := totalexc + object$resultmean[, get(paste0(object$str_exc[i], "_k"))] * get(object$str_exc[i])]
      }
      newdata[, paste0(object$str_signal, "_estimated") := generate.1order(time = get(object$str_time),
                                                                           excitation = totalexc,
                                                                           y0 = object$resultmean[, yeq],
                                                                           tau = object$resultmean[, tau],
                                                                           k = 1,
                                                                           yeq = object$resultmean[, yeq])$y]
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
  res$str_time <- object$str_time
  res$str_exc <- object$str_exc
  res$str_signal <- object$str_signal
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
