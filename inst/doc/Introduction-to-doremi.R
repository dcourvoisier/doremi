## ------------------------------------------------------------------------
devtools::load_all("C:/Users/Alba/Documents/Unige/doremi")

## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  warnings = FALSE,
  comment = "#>"
)

## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
mydataa1 <- simulation_generate_order1(nindividuals = 4, 
                                    dampingtime = 10, 
                                    amplitude = 1, 
                                    nexc = 3, 
                                    duration = 10, 
                                    deltatf = 0.5,
                                    tmax = 100,
                                    minspacing = 20,
                                    internoise = 0, 
                                    intranoise = 0)

## ------------------------------------------------------------------------
head(mydataa1)

## ----fig.width = 7,fig.height = 6----------------------------------------
ggplot2::ggplot( data = mydataa1$data ) +
  ggplot2::geom_point(ggplot2::aes(timecol,dampedsignalraw, colour = "Signal-no noise"))+
  ggplot2::geom_point(ggplot2::aes(timecol,dampedsignal, colour = "Signal with 0% intra-noise"))+
  ggplot2::geom_line(ggplot2::aes(timecol,excitation,colour = "Excitation"))+
  ggplot2::facet_wrap(~id)+
  ggplot2::labs(x = "Time (unit)",
           y = "Signal (unit)",
           colour = "")

## ------------------------------------------------------------------------
# Generation of signals with intra and inter-noise
mydataa2 <- simulation_generate_order1(nindividuals = 4, 
                                    dampingtime = 10, 
                                    amplitude = 1, 
                                    nexc = 3, 
                                    duration = 10, 
                                    deltatf = 0.5,
                                    tmax = 100,
                                    minspacing = 20,
                                    internoise = 0.4, 
                                    intranoise = 0.2)

## ----fig.width = 7,fig.height = 6----------------------------------------
ggplot2::ggplot( data = mydataa2$data ) +
  ggplot2::geom_point(ggplot2::aes(timecol,dampedsignalraw, colour = "Signal-no noise"))+
  ggplot2::geom_point(ggplot2::aes(timecol,dampedsignal, colour = "Signal with 20% intra-noise"))+
  ggplot2::geom_line(ggplot2::aes(timecol,excitation,colour = "Excitation"))+
  ggplot2::facet_wrap(~id)+
  ggplot2::labs(x = "Time (unit)",
           y = "Signal (unit)",
           colour = "")

## ------------------------------------------------------------------------
resultb1 <- doremi_analyse_order1(userdata = mydataa2$data,
                                id = "id",
                                input ="excitation",
                                time ="timecol",
                                signal = "dampedsignal",
                                embedding = 5)

## ------------------------------------------------------------------------
head(resultb1$data)

## ------------------------------------------------------------------------
resultb1$resultID

## ------------------------------------------------------------------------
resultb1$resultmean

## ------------------------------------------------------------------------
resultb1$regression

## ------------------------------------------------------------------------
resultb1$estimated

## ----fig.width = 7, fig.height = 6---------------------------------------
plot(resultb1)

## ----fig.width = 7, fig.height = 6---------------------------------------
plot(resultb1, id = 3)
plot(resultb1, id = c(1,5))

## ------------------------------------------------------------------------
#Simulating data with these hypothesis
#Generating the three excitation signals:
e1 <- excitation_function ( amplitude = 10, 
                            nexc = 1, 
                            duration = 10, 
                            deltatf = 1, 
                            tmax = 100,
                            minspacing = 20)
e2 <- excitation_function ( amplitude = 10, 
                            nexc = 1, 
                            duration = 10, 
                            deltatf = 1, 
                            tmax = 100,
                            minspacing = 20)
e3 <- excitation_function ( amplitude = 10, 
                            nexc = 1, 
                            duration = 10, 
                            deltatf = 1, 
                            tmax = 100,
                            minspacing = 20)
et1 <- e1$y+3*e2$y+5*e3$y
timt1 <- e3$t  #we can use any of the three time vectors as they are identical for the three excitations
y1 <- doremi_generate_order1(10,et1,timt1)$y

#Signals for the second individual
et2 <- e1$y+2.5*e2$y+4*e3$y
y2 <- doremi_generate_order1(10,et2,timt1)$y

#Generating table with signals
mydatab2 <- data.table::setDT(list(id = rep(c(1, 2), c(length(et1), length(et2))), 
                       timecol = c(timt1, timt1),
                       excitation1 = c(e1$y, e1$y),
                       excitation2= c(e2$y, e2$y),
                       excitation3 = c(e3$y, e3$y),
                       excitation = c(et1, et2), 
                       signalcol = c(y1, y2)))

## ----fig.width = 7,fig.height = 4----------------------------------------
#Plotting signals
ggplot2::ggplot( data = mydatab2 ) +
  ggplot2::geom_line(ggplot2::aes(timecol,signalcol, colour = "Signal-no noise"))+
  ggplot2::geom_line(ggplot2::aes(timecol,excitation,colour = "Excitation"))+
  ggplot2::facet_wrap(~id)+
  ggplot2::labs(x = "Time (unit)",
           y = "Signal (unit)",
           colour = "")

## ------------------------------------------------------------------------
#Analysing signals
resultb2 <- doremi_analyse_order1(userdata = mydatab2,
                                id = "id",
                                input = c("excitation1", "excitation2", "excitation3"),
                                time ="timecol",
                                signal = "signalcol",
                                embedding = 5)

#Looking for the calculation of the coefficients of the excitation
resultb2$resultid

## ----fig.width = 7,fig.height = 4----------------------------------------
#Plotting signals
plot(resultb2)

## ------------------------------------------------------------------------
#Simulating data with these hypothesis
mydatab3 <- simulation_generate_order1(nindividuals = 6, 
                                      dampingtime = 10, 
                                      amplitude = -1,
                                      nexc = 1, 
                                      duration = 50, 
                                      deltatf = 1,
                                      tmax = 50,
                                      minspacing = 0, 
                                      internoise = 0.4, 
                                      intranoise = 0.2)
mydatab3$data[, dampedsignal := dampedsignal + 80]
mydatab3$data[, dampedsignalraw := dampedsignalraw + 80] #To ensure that there are no negative values when fitting the log


## ------------------------------------------------------------------------
#Analysing
resultb3 <- doremi_analyse_order1(userdata = mydatab3$data,
                                id = "id",
                                time = "timecol",
                                signal = "dampedsignal",
                                embedding = 5)

## ----fig.width = 7, fig.height = 6---------------------------------------
#Plotting
plot(resultb3)

## ------------------------------------------------------------------------
#Creating the data table
mydatab4 <- data.table::setDT(list(timecol = timt1,
                       excitation = et1, 
                       signalcol = y1))

## ------------------------------------------------------------------------
#Analysing
resultb4 <- doremi_analyse_order1(userdata = mydatab4,
                                input = "excitation",
                                time ="timecol",
                                signal = "signalcol",
                                embedding = 5)

## ----fig.width = 5,fig.height = 4, fig.pos = 0.5-------------------------
#Plotting 
plot(resultb4)

## ------------------------------------------------------------------------
mydatab5 <- simulation_generate_order1(nindividuals = 4, 
                                    dampingtime = 10, 
                                    amplitude = 1, 
                                    nexc = 3, 
                                    duration = 10, 
                                    deltatf = 0.5,
                                    tmax = 100,
                                    minspacing = 20,
                                    internoise = 0.1, 
                                    intranoise = 0.2)
#Keeping half of the rows selected randomly
mydatab5$data[, keep := sample(c(0,1), .N, replace = TRUE), by = id]
mydatab5rd <- mydatab5$data[keep == 1]

#Analysing the resulting signal
resultb5 <- doremi_analyse_order1(userdata = mydatab5rd,
                                input = "excitation",
                                time ="timecol",
                                signal = "dampedsignal",
                                embedding = 5)

## ----fig.width = 7,fig.height = 6----------------------------------------
ggplot2::ggplot( data = mydatab5$rawdata ) +
  ggplot2::geom_point(ggplot2::aes(timecol,dampedsignalraw, colour = "Pseudo-continuous signal. No noise"))+
  ggplot2::geom_line(ggplot2::aes(timecol,excitation,colour = "Excitation"))+
  ggplot2::geom_point(data = mydatab5rd, ggplot2::aes(timecol,dampedsignalraw, colour = "Random sampling"))+
  ggplot2::facet_wrap(~id)+
  ggplot2::labs(x = "Time (unit)",
           y = "Signal (unit)",
           colour = "")

## ----fig.width = 7, fig.height = 6---------------------------------------
#plot(resultb5)
resultb5$data

## ------------------------------------------------------------------------
resultc1 <- doremi_analyse_order1(userdata = cardio,
                                id = "id",
                                input = "load",
                                time ="timecol",
                                signal = "hr",
                                embedding = 5)

## ----fig.width = 10,fig.height = 10--------------------------------------
plot(resultc1, id = 1:21)

## ------------------------------------------------------------------------
resultc2 <- doremi_analyse_order1(userdata = rotation,
                                id = "id",
                                time ="numtrials",
                                signal = "meanRT",
                                embedding = 5)

## ----fig.width = 10, fig.height = 10-------------------------------------
plot(resultc2, id = 1:17)

## ----fig.width = 5,fig.height = 4, fig.pos = 0.5-------------------------
#Input data
e1 <- excitation_function ( amplitude = 10,
                            nexc = 1,
                            duration = 10,
                            deltatf = 1,
                            tmax = 100,
                            minspacing = 20)
e2 <- excitation_function ( amplitude = 10,
                            nexc = 1,
                            duration = 10,
                            deltatf = 1,
                            tmax = 100,
                            minspacing = 20)
e3 <- excitation_function ( amplitude = 10,
                            nexc = 1,
                            duration = 10,
                            deltatf = 1,
                            tmax = 100,
                            minspacing = 20)
et1 <- e1$y+3*e2$y+5*e3$y
timt1 <- e3$t  #we can use any of the three time vectors as they are identical for the three excitations
y1 <- doremi_generate_order1(10,et1,timt1)$y

#Signals for the second individual
et2 <- e1$y+2.5*e2$y+4*e3$y
y2 <- doremi_generate_order1(10,et2,timt1)$y

#Create data table with the pair e1,y1
mydata<-setDT(list(tim=timt1,exc=et1,y1=y1))

#Analysis 
resultd1 <- doremi_analyse_order1(mydata,
                                 input = "exc",
                                 time = "tim",
                                 signal = "y1",
                                 embedding = 5)
                                 
#Create data table with et2 that will be supplied as new excitation for the predict function
myexc_table <- mydata[, exc:=et2]

#Calling the predict function
predresultd1<- predict (resultd1,
                       newdata = myexc_table)

#Compare with calculated signal y2
q<-plot(predresultd1)
q+geom_point(aes(predresultd1$data$tim, y2,colour = "y2"))

