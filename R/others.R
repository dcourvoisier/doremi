# Call to global variables ------------------------------------------------

globalVariables(c("inputcol", "timecol", "excitation_sum", "expfit_A"))
globalVariables(c("dampedsignalraw", "amplitudenorm", "dampedsignal", "exc_max", "exc_min", "excitation", "id", "intervals", "totalexc", "ymax", "ymin"))

#Classes and methods
#Classes documented directly in the analysis function (object "doremi")
#'@import ggplot2
plot.doremi=function(doremi_object,
                     id = 1:5,
                     signalcolumn = "dampedsignal",
                     timecolumn = "timecol",
                     exc = NULL,
                     raw = NULL,
                     xlabel = "Time",
                     ylabel = "Signal"){
  if (is.null(doremi_object$data$id)){ #Single individual. Doremi object doesn't have id column
    facets<-as.character(id)
    p <- ggplot(doremi_object$data[doremi_object$data$id %in% facets])
  }else{# Multiple individuals
    p <- ggplot(doremi_object$data)
  }
  P <- p+geom_point(aes(get(timecolumn), get(signalcolumn), colour = "Signal")) +
  facet_wrap(~id,scales = "free")+
  geom_point(data = doremi_object$estimated[doremi_object$estimated$id %in% facets], aes(timecol, ymin, colour = "Signal estimated - ymin")) +
  geom_point(estimated = doremi_object$estimated[doremi_object$estimated$id %in% facets], aes(timecol, ymax, colour = "Signal estimated - ymax"))
  if(!is.null(exc)){
    p<-p+geom_line(aes(get(timecolumn), get(exc), colour = "Excitation"))
  }
  if(!is.null(raw)){
    p<-p+geom_line(aes(get(timecolumn), get(raw), colour= "Signal- no noise"))
  }
  p+labs(x = xlabel,
         y = ylabel,
         colour = "Legend")
}
