## ----options, include = FALSE-------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----libraries, setup---------------------------------------------------------
library(doremi)
library(ggplot2)
library(data.table)
library(magrittr)
library(deSolve)
library(gridExtra)
set.seed(1)

## ----different forms according to xi,fig.width = 5, fig.height = 4, fig.align = "center"----
data11a <- rbindlist(lapply(seq(0,2,0.2), function(eps){
  generate.2order(time = 0:49, y0 = 1, xi = eps, period = 30)[,xi := eps][]
}))
# plot
ggplot2::ggplot(data11a,ggplot2::aes(t,y,color = as.factor(xi)))+
  ggplot2::geom_line() +
  ggplot2::labs(x = "time (arb. unit)", y = "signal (arb. unit)", colour = "xi")

## ----simulation example 1-----------------------------------------------------
time <- 0:100
data1 <- generate.panel.2order(time = time,
                               y0 = 10,
                               v0 = 0,
                               xi = 0.1,
                               period = 30,
                               k = 1,
                               yeq = 2,
                               nind = 6,
                               internoise = 0.3,
                               intranoise = 5)
data1

## ----dlo plot data1,fig.width = 7, fig.height = 6, fig.align = "center"-------
plot(data1) +
  ggplot2::geom_hline(yintercept=0)

## ----simulation example2------------------------------------------------------
data2a <- generate.panel.2order(time = 0:99,
                               y0 = 1,
                               v0 = 0,
                               xi = 0,
                               period = 30,
                               k = 1,
                               yeq = 0,
                               nind = 1,
                               internoise = 0,
                               intranoise = 5)
data2b <- generate.panel.2order(time = 0:99,
                               y0 = 1,
                               v0 = 0,
                               xi = 1,
                               period = 30,
                               k = 1,
                               yeq = 0,
                               nind = 1,
                               internoise = 0,
                               intranoise = 5)
data2c <- generate.panel.2order(time = 0:99,
                               y0 = 1,
                               v0 = 0,
                               xi = 2,
                               period = 30,
                               k = 1,
                              yeq = 0,
                               nind = 1,
                               internoise = 0,
                               intranoise = 5)


## ----plot example2,fig.width = 7, fig.height = 5, fig.align = "center"--------
grid.arrange(plot(data2a)+ggplot2::ggtitle("undamped, xi=0"), plot(data2b)+ggplot2::ggtitle("critically damped, xi=1"), plot(data2c)+ggplot2::ggtitle("overdamped, xi=2"), ncol= 3)

## ----analysis, example1-------------------------------------------------------
res1 <- analyze.2order(data = data1,
                        id = "id",
                        time ="time",
                        signal = "signal",
                        dermethod = "glla",
                        derparam = 9,
                        order = 4)

## ----analysis plot res1,fig.width = 7, fig.height = 6, fig.align = "center"----
plot(res1)

## ----print res1 and components------------------------------------------------
res1$resultid

## ----print res1---------------------------------------------------------------
res1

## ----2nd order with excitation example1---------------------------------------
time <- 0:100
data1 <- generate.panel.2order(time = time,
                               excitation = as.numeric(time>20),
                               y0 = 0,
                               v0 = 0,
                               xi = 0.1,
                               period = 30,
                               k = 1,
                               yeq = 0,
                               nind = 5,
                               internoise = 0,
                               intranoise = 0)

## ----plot data1,fig.width = 7, fig.height = 6, fig.align = "center"-----------
plot(data1)

## ----example2-----------------------------------------------------------------
# Generation of signals with intra and inter-noise
time <- 0:100
data2 <- generate.panel.2order(time = time,
                               excitation = as.numeric(time>20),
                               y0 = 0,
                               v0 = 0,
                               xi = 0.1,
                               period = 30,
                               k = 1,
                               yeq = 0,
                               nind = 5,
                               internoise = 0.3,
                               intranoise = 5)

## ----plot data2,fig.width = 7, fig.height = 6, fig.align = "center"-----------
plot(data2)

## ----example3-----------------------------------------------------------------
time <- 0:99
data3 <- generate.2order(time = time,
                               excitation = as.numeric(time>20),
                               y0 = 0,
                               v0 = 0,
                               t0 = 0,
                               exc0 = 0,
                               xi = 0.1,
                               period = 30,
                               k = 1,
                               yeq = 0)
data3[t==43]

## ----example3 rebuilding signal-----------------------------------------------
data3_t0_43 <- generate.2order(time = time,
                               excitation = as.numeric(time>20),
                               y0 = 1.077549,
                               v0 = -1.311226e-01,
                               t0 = 43,
                               exc0 = 1,
                               xi = 0.1,
                               period = 30,
                               k = 1)


## ----plot example3,fig.width = 6, fig.height = 4, fig.pos = 0.5, fig.align = "center"----
ggplot2::ggplot() + 
  ggplot2::geom_line(data=data3,ggplot2::aes(t,y,color="Solution with t0=0")) +
  ggplot2::geom_point(data=data3_t0_43,ggplot2::aes(t,y,color="Solution with t0=43"))

## ----simulation example4------------------------------------------------------
t <-0:99
xi <- 0.2
period <- 30
dlo <- generate.2order(time = t,y0 = 1, xi = xi, period = period)
driven_dlo <- generate.2order(time = t, excitation = 5*sin(t), y0 = 1, xi = xi, period = period)

## ----simulation plot example4, fig.width = 5, fig.height = 4, fig.align = "center"----
ggplot2::ggplot() + 
  ggplot2::geom_point(data=dlo,ggplot2::aes(t,y,color="dlo")) + 
  ggplot2::geom_line(data=driven_dlo,ggplot2::aes(t,y,color="driven dlo"))+
  labs(colour = "") +
  ggtitle("Dlo versus driven dlo (excitation=5sin(t))")+
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5))

## ----analysis example1--------------------------------------------------------
res1 <- analyze.2order(data = data1[id==1],
                      input = "excitation",
                      time ="time",
                      signal = "signal",
                      dermethod = "fda",
                      derparam = 0.8,
                      verbose=T)

## ----plot res1,fig.width = 6, fig.height = 4, fig.pos = 0.5, fig.align = "center"----
plot(res1)

## ----res2---------------------------------------------------------------------
res2 <- analyze.2order(data = data2,
                        id = "id",
                        input ="excitation",
                        time ="time",
                        signal = "signal",
                        dermethod = "gold",
                        derparam = 5,
                        order = 4)

## ----plot res2,fig.width = 7, fig.height = 6, fig.align = "center"------------
plot(res2)

## ----res3---------------------------------------------------------------------
res3<-optimum_param (data=data2,
                      id="id",
                      input="excitation",
                      time="time",
                      signal="signal",
                      model = "2order",
                      dermethod = "glla",
                      order = 2,
                      pmin = 5,
                      pmax = 17,
                      pstep = 2)
res3$analysis
res3$summary_opt
res3$d

## ----plot res3, fig.width = 7, fig.height = 6, fig.align = "center"-----------
plot(res3)

## ----optimum analysis res3b---------------------------------------------------
res3b <- analyze.2order(data = data2,
                        id = "id",
                        input ="excitation",
                        time ="time",
                        signal = "signal",
                        dermethod = "glla",
                        derparam = 13,
                        order = 2)
res3b

## ----plot res3b,fig.width = 7, fig.height = 6, fig.align = "center"-----------
plot(res3b)

## ----example4-----------------------------------------------------------------
#Simulating data with these hypothesis
#Generating the three excitation signals:
u1 <- generate.excitation (amplitude = 10, 
                           nexc = 1, 
                           duration = 10, 
                           deltatf = 1, 
                           tmax = 100,
                           minspacing = 20)
u2 <- generate.excitation (amplitude = 10, 
                           nexc = 1, 
                           duration = 10, 
                           deltatf = 1, 
                           tmax = 100,
                           minspacing = 20)
u3 <- generate.excitation (amplitude = 10, 
                           nexc = 1, 
                           duration = 10, 
                           deltatf = 1, 
                            tmax = 100,
                           minspacing = 20)
# Arbitrarily choosing a = 1, b = 2 and c = 5 for the first individual
et1 <- u1$exc + 3 * u2$exc + 5 * u3$exc
timt1 <- u3$t  #we can use any of the three time vectors as they are identical for the three excitations
y1 <- generate.2order(time = timt1,
                            excitation = et1)$y
#as we are using the $y argument of the object generated

#Signals for the second individual;
# Arbitrarily choosing a = 1, b = 2.5 and c = 4 for the second individual
et2 <- u1$exc + 2.5 * u2$exc + 4 * u3$exc
y2 <- generate.2order(time = timt1,
                            excitation = et2)$y 

#Generating table with signals
dataa4 <- data.frame(id = rep(c(1, 2), c(length(et1), length(et2))), 
                 time = c(timt1, timt1),
                 excitation1 = c(u1$exc, u1$exc),
                 excitation2 = c(u2$exc, u2$exc),
                 excitation3 = c(u3$exc, u3$exc),
                 excitation = c(et1, et2), 
                 signalcol = c(y1, y2))

## ----plot example4,fig.width = 7, fig.height = 4, fig.align = "center"--------
#Plotting signals
ggplot2::ggplot( data = dataa4) +
  ggplot2::geom_line(ggplot2::aes(time,signalcol, colour = "Signal-no noise"))+
  ggplot2::geom_line(ggplot2::aes(time,excitation,colour = "Total excitation"))+
  ggplot2::facet_wrap(~id)+
  ggplot2::labs(x = "Time (s)",
           y = "Signal (arb. unit)",
           colour = "")

## ----res4---------------------------------------------------------------------
#Analyzing signals
res4 <- analyze.2order(data = dataa4,
                       id = "id",
                       input = c("excitation1", "excitation2", "excitation3"),
                       time = "time",
                       signal = "signalcol",
                       dermethod = "fda",
                       derparam = 0.1)

#Looking for the calculation of the coefficients of the excitation
res4
res4$resultid


## ----plot res4,fig.width = 7, fig.height = 4, fig.align = "center"------------
#Plotting signals
plot(res4)

