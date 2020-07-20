## ----options, include=FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = FALSE,
  comment = "#>"
)

## ----libraries, setup---------------------------------------------------------
library(doremi)
library(ggplot2)
library(data.table)
library(futile.logger)
set.seed(1)

## ----simulation example1a-----------------------------------------------------
res1a <- generate.panel.1order(time = 0:99,
                               excitation = NULL,
                               y0 = 1,
                               t0 = 0,
                               exc0 = 0,
                               tau = 10,
                               k = 1,
                               yeq = 0,
                               nind = 4,
                               internoise = 0,
                               intranoise = 0)

## ----res1a--------------------------------------------------------------------
res1a

## ----plot res_1a, fig.align="center", fig.height=6, fig.width=7---------------
plot(res1a)

## ----changing initial condition res1b-----------------------------------------
res1b <- generate.panel.1order(time = 0:49,
                               excitation = NULL,
                               y0 = 3,
                               t0 = 0,
                               exc0 = 0,
                               tau = 5,
                               k = 0,
                               yeq = 1.5,
                               nind = 4,
                               internoise = 0,
                               intranoise = 0)

## ----changing initial condition plot res1b,fig.width = 7, fig.height = 6, fig.align = "center"----
plot(res1b) + scale_y_continuous(limits = c(0, 3))

## ----noise res2a--------------------------------------------------------------
# Generation of signals with intra and inter-noise
res2a <- generate.panel.1order(time = 0:49,
                               excitation = NULL,
                               y0 = 1,
                               t0 = 0,
                               exc0 = 0,
                               tau = 5,
                               k = 1,
                               yeq = 0,
                               nind = 6,
                               internoise = 0.4,
                               intranoise = 5)

## ----noise plot res2a, fig.width = 7, fig.height = 6, fig.align = "center"----
plot(res2a)

## ----analysis example3--------------------------------------------------------
#Simulating data with these hypothesis
data3 <- generate.panel.1order(time = 0:50,
                               excitation = NULL,
                               y0 = 0.5,
                               t0 = 0,
                               exc0 = 0,
                               tau = 10,
                               k = 1,
                               yeq = 0,
                               nind = 3,
                               internoise = 0.2,
                               intranoise = 100)

## ----analysis res3------------------------------------------------------------
#Analyzing

res3 <- analyze.1order(data = data3,
                      id = "id",
                      time = "time",
                      signal = "signal",
                      dermethod = "gold",
                      derparam = 11)

## ----analysis plot res3,fig.width = 7, fig.height = 6, fig.align = "center"----
#Plotting
plot(res3)

## ----excitation term example1a------------------------------------------------
U1a <- generate.excitation(amplitude = 1,
                           nexc = 3,
                           duration = 10,
                           deltatf = 1,
                           tmax = 100,
                           minspacing = 20)
                               
res1a <- generate.panel.1order(time = U1a$t,
                               excitation = U1a$exc,
                               y0 = 0,
                               t0 = 0,
                               exc0 = 0,
                               tau = 10,
                               k = 1,
                               yeq = 0,
                               nind = 4,
                               internoise = 0,
                               intranoise = 0)

## ----changing initial condition-----------------------------------------------
U1b <- generate.excitation(amplitude = 1,
                           nexc = 1,
                           duration = 50,
                           deltatf = 1,
                           tmax = 100,
                           minspacing = 20)
                               
res1b <- generate.panel.1order(time = U1b$t,
                               excitation = U1b$exc,
                               y0 = 3,
                               t0 = 0,
                               exc0 = 0,
                               tau = 10,
                               k = 5,
                               yeq = 2,
                               nind = 4,
                               internoise = 0,
                               intranoise = 0)

## ----plot res1b,fig.width = 7, fig.height = 6, fig.align = "center"-----------
plot(res1b)

## ----example2-----------------------------------------------------------------
# Generation of signals with intra and inter-noise
U2a <- generate.excitation(amplitude = 1,
                           nexc = 3,
                           duration = 10,
                           deltatf = 1,
                           tmax = 100,
                           minspacing = 20)
res2a <- generate.panel.1order(time = U2a$t,
                               excitation = U2a$exc,
                               y0 = 0,
                               t0 = 0,
                               exc0 = 0,
                               tau = 10,
                               k = 1,
                               yeq = 0,
                               nind = 4,
                               internoise = 0.4,
                               intranoise = 5)

## ----plot res2a,fig.width = 7, fig.height = 6, fig.align = "center"-----------
plot(res2a)

## ----example3-----------------------------------------------------------------
data3 <- data.frame(U2a, 
                 signalcol = res2a[id==1]$signal)

## ----res3---------------------------------------------------------------------
#Analyzing
res3 <- analyze.1order(data = data3,
                      input = "exc",
                      time ="t",
                      signal = "signalcol",
                      verbose=T)

## ----plot res3,fig.width = 6, fig.height = 4, fig.pos = 0.5, fig.align = "center"----
#Plotting 
plot(res3)

## ----analysis res3a-----------------------------------------------------------
res3a <- analyze.1order(data = res2a,
                        id = "id",
                        input ="excitation",
                        time ="time",
                        signal = "signal",
                        dermethod = "gold",
                        derparam = 5,
                        order = 1)

## ----analysis res3a print-----------------------------------------------------
res3a 

## ----summary res3a------------------------------------------------------------
summary(res3a) 

## ----head res3a---------------------------------------------------------------
head(res3a$data)

## ----components of res3a------------------------------------------------------
res3a$regression

## ----resultmean---------------------------------------------------------------
res3a$resultmean

## ----resultid-----------------------------------------------------------------
res3a$resultid

## ----plot res3a,fig.width = 7, fig.height = 6, fig.align = "center"-----------
plot(res3a)

## ----more plot res3a,fig.width = 7, fig.height = 6, fig.align = "center"------
plot(res3a, id = 3)
plot(res3a, id = c(1,4))

## ----example 3b---------------------------------------------------------------
time <- 0:100
simu_data <- generate.panel.1order(time = time,
                                   excitation = as.numeric(time>40),
                                   y0 = 0,
                                   t0 = 0,
                                   exc0 = 0,
                                   tau = 10,
                                   k = 1,
                                   yeq = 0,
                                   nind = 10,
                                   internoise = 0.2,
                                   intranoise = 5)

## ----plot,fig.width = 7, fig.height = 6, fig.align = "center"-----------------
plot(simu_data)

## ----res3b--------------------------------------------------------------------
res3b<-optimum_param (data = as.data.table(simu_data),
                      id = "id",
                      input = "excitation",
                      time = "time",
                      signal = "signal",
                      model = "1order",
                      dermethod = "glla",
                      pmin = 3,
                      pmax = 11,
                      pstep = 2)
res3b$analysis
res3b$summary_opt
res3b$d

## ----plot res3b,fig.width = 7, fig.height = 6, fig.align = "center"-----------
plot(res3b)

## ----example4-----------------------------------------------------------------
#Simulating data with these hypothesis
#Generating the three excitation signals:
t <- 0:99
u1 <- as.numeric(t>20 & t<40)
u2 <- as.numeric(t>50 & t<60)
u3 <- as.numeric(t>80 & t<85)
# Arbitrarily choosing a = 1, b = 2 and c = 5 for the first individual
et1 <- u1 + 3 * u2 + 5 * u3
y1 <- generate.1order(time = t,
                      excitation = et1,
                      y0 = 0,
                      t0 = 0,
                      exc0 = 0,
                      tau = 10,
                      k = 1,
                      yeq = 0)$y
#as we are using the $y argument of the object generated

#Signals for the second individual;
# Arbitrarily choosing a = 1, b = 2.5 and c = 4 for the second individual
et2 <- u1 + 2.5 * u2 + 4 * u3
y2 <- generate.1order(time = t,
                      excitation = et2,
                      y0 = 0,
                      t0 = 0,
                      exc0 = 0,
                      tau = 10,
                      k = 1,
                      yeq = 0)$y 

#Generating a table with the signals
data4 <- data.frame(id = rep(c(1, 2), c(length(et1), length(et2))), 
                 time = c(t, t),
                 excitation1 = c(u1, u1),
                 excitation2 = c(u2, u2),
                 excitation3 = c(u3, u3),
                 excitation = c(et1, et2), 
                 signalcol = c(y1, y2))

## ----plot example4,fig.width = 7, fig.height = 4, fig.align = "center"--------
#Plotting signals
ggplot2::ggplot( data = data4) +
  ggplot2::geom_line(ggplot2::aes(time,signalcol, colour = "Signal-no noise"))+
  ggplot2::geom_line(ggplot2::aes(time,excitation,colour = "Total excitation"))+
  ggplot2::facet_wrap(~id)+
  ggplot2::labs(x = "Time (s)",
           y = "Signal (arb. unit)",
           colour = "")

## ----res4---------------------------------------------------------------------
#Analyzing signals
res4 <- analyze.1order(data = data4,
                       id = "id",
                       input = c("excitation1", "excitation2", "excitation3"),
                       time = "time",
                       signal = "signalcol",
                       dermethod = "fda",
                       derparam = 0.1)

#Looking for the calculation of the coefficients of the excitation
res4
res4$resultid


## ----fig.width = 7, fig.height = 4, fig.align = "center"----------------------
#Plotting signals
plot(res4)

## ----example5-----------------------------------------------------------------
t <- 0:200
data5 <- generate.panel.1order(time = t,
                               excitation = as.numeric(t>50 & t<100),
                               y0 = 0,
                               t0 = 0,
                               exc0 = 0,
                               tau = 10,
                               k = 1,
                               yeq = 0,
                               nind = 6,
                               internoise = 0.4,
                               intranoise = 100)

## ----plot data5, fig.width = 7, fig.height = 6, fig.align = "center"----------
plot(data5)

## ----missing data-------------------------------------------------------------
#Keeping one third of the rows selected randomly from the full data set
data5rd <- as.data.table(data5[sample(nrow(data5), nrow(data5)/3), ])
data5rd <- data5rd[order(id,time)]

## ----plot missing data,fig.width = 7, fig.height = 6, fig.align = "center"----
ggplot2::ggplot( data = data5 ) +
  ggplot2::geom_point(ggplot2::aes(time,signalraw, colour = "Full data set. No noise"))+
  ggplot2::geom_line(ggplot2::aes(time,excitation,colour = "Excitation"))+
  ggplot2::geom_point(data = data5rd, ggplot2::aes(time,signalraw, colour = "Random sampling"))+
  ggplot2::facet_wrap(~id)+
  ggplot2::labs(x = "Time (unit)",
           y = "Signal (unit)",
           colour = "")

## ----res7---------------------------------------------------------------------
res7 <- analyze.1order(data = data5rd,
                       id = "id",
                       input = "excitation",
                       time ="time",
                       signal = "signal")

## ----plot res7,fig.width = 7, fig.height = 6, fig.align = "center"------------
ggplot2::ggplot( data = data5 ) +
  ggplot2::geom_line(ggplot2::aes(time,signalraw, colour = "Original signal"))+
  ggplot2::geom_line(ggplot2::aes(time,excitation,colour = "Excitation"))+
  ggplot2::geom_point(data = data5rd, ggplot2::aes(time,signal, colour = "Random sampled signal"))+
  ggplot2::geom_line(data = res7$data, ggplot2::aes(time,signal_estimated, colour = "Estimated signal. Missing data points"))+
  ggplot2::facet_wrap(~id)+
  ggplot2::labs(x = "Time (unit)",
           y = "Signal (unit)",
           colour = "") +
  ggplot2::theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5))



## ----cardio data--------------------------------------------------------------
resc1a <- analyze.1order(data = cardio,
                 id = "id",
                 input = "load",
                 time ="time",
                 signal = "hr",
                 dermethod = "gold",
                 derparam = 5)

## ----plot cardio data,fig.width = 7, fig.height = 16, fig.align = "center"----
plot(resc1a, id = 1:21) + ggplot2::facet_wrap(~id,ncol=3,scales="free")

## ----cardio multiple excitations----------------------------------------------
mydata<-cardio

#Decomposing the load
#Maximum number of extra columns to create
Ncol <- max(mydata[,length(unique(load[load != 0])),by = id]$V1)
emptycols<-c(paste0("input",1:Ncol))

# Creating column that contains unique values in load greater than 0 for decomposition
mydata[,Ncol := length(unique(load[load != 0])),by = id]

#Filling the created empty columns with the values
for (i in 1:length(unique(mydata$id))){
  Ncol<-mydata[id==unique(mydata$id)[i],length(unique(load[load != 0]))]
  mydata[id==unique(mydata$id)[i],emptycols[1:Ncol] := lapply(unique(load[load != 0]),function(x){ ifelse(load == x,x,0)}), by = id]
}

resc1b <- analyze.1order(data = mydata,
                        id = "id",
                        input = emptycols,
                        time = "time",
                        signal = "hr",
                        dermethod = "gold",
                        derparam = 5)

## ----plot cardio mult excitations,fig.width = 7, fig.height = 16, fig.align = "center"----
plot(resc1b,id=1:21)+ ggplot2::facet_wrap(~id,ncol=3,scales="free")

## ----mental rotation data-----------------------------------------------------
dermethod<- "fda"
pmin = 0.1
pmax = 1
pstep = 0.1

restemp <- optimum_param(data = rotation,
                          id = "id",
                          time ="days",
                          signal = "meanRT",
                          dermethod = dermethod,
                          model = "1order",
                          pmin = pmin,
                          pmax = pmax,
                          pstep = pstep)
restemp$summary_opt
resc2a <- analyze.1order(data = rotation,
                 id = "id",
                 time ="days",
                 signal = "meanRT",
                 dermethod = dermethod,
                 derparam = restemp$d)

## ----plot menta rotation data,fig.width = 7, fig.height = 6, fig.align = "center"----
plot(resc2a, id = 1:17)

## ----print summary predict, fig.width = 5,fig.height = 4, fig.pos = 0.5, fig.align = "center"----
#Input data from previous example
datad1<-setDT(data4)[id==1,.(id,time,excitation,signalcol)]

#Analysis 
resd1 <- analyze.1order(data = datad1,
                        input = "excitation",
                        time = "time",
                        signal = "signalcol",
                        dermethod = "glla",
                        derparam = 7)
                                 
#Creating data frame with et2 that will be supplied as new excitation for the predict function
datad1n <- copy(datad1)
datad1n$excitation <- as.numeric(datad1$time>50 & datad1$time<60)


#Calling the predict function
predresd1<- predict(resd1,
                    newdata = datad1n)

#Comparing predicted value with signal y2
ggplot2::ggplot(data = predresd1$data) +
ggplot2::geom_point(ggplot2::aes(time, signalcol_estimated,colour = "predicted signalcol")) +
ggplot2::geom_point(ggplot2::aes(time, signalcol,colour = "signalcol")) +
ggplot2::geom_line(data= datad1, ggplot2::aes(time, excitation,colour = "Original excitation")) + 
ggplot2::geom_line(data= datad1n, ggplot2::aes(time, excitation,colour = "New excitation"))  


