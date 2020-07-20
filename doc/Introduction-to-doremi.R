## ---- include = FALSE---------------------------------------------------------
 library(ggplot2)
 library(data.table)
 library(doremi)
set.seed(1)
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = FALSE,
  comment = "#>"
)

## ----fig.width = 5, fig.height = 4, fig.align = "center", echo = FALSE--------
time <- 0:90
exc <- rep(c(0,1,0),c(11,30,50))
signal<- generate.1order(time = time,
                         excitation = exc,
                         y0 = 0,
                         t0 = 0,
                         exc0 = 0,
                         tau = 5,
                         k = 1,
                         yeq = 0)
signal<-cbind(data.table::as.data.table(signal),exc)
ggplot2::ggplot(data = as.data.table(signal)) +
ggplot2::ggtitle( "First order differential equation solution")+
  ggplot2::geom_line(ggplot2::aes(t,y, colour = "Signal"))+
  ggplot2::geom_line(ggplot2::aes(t,exc, colour = "Excitation"))+
  ggplot2::geom_hline(yintercept=0.63*max(signal$y), linetype="dashed", colour = "gray")+
  ggplot2::geom_hline(yintercept=0.37*max(signal$y), linetype="dashed", colour = "gray")+
  ggplot2::geom_vline(xintercept=signal$t[signal$y==max(signal$y)], colour = "gray")+
  ggplot2::geom_vline(xintercept=50, colour = "gray")+
  ggplot2::geom_vline(xintercept=19, colour = "gray")+
  ggplot2::geom_vline(xintercept=10, colour = "gray")+
  ggplot2::annotate("segment", x = 10, xend = 19, y = -0.1, yend = -0.1, colour = "dark green", size = 1)+
  ggplot2::annotate("text", x = 15, y = -0.2, label = "tau", parse = TRUE, colour = "dark green")+
  
  ggplot2::annotate("segment", x = 40, xend = 50, y = -0.1, yend = -0.1, colour = "dark green", size = 1)+
  ggplot2::annotate("text", x = 45, y = -0.2, label = "tau", parse = TRUE, colour = "dark green")+

  ggplot2::annotate("text", x = 75, y = 0.63*max(signal$y), label = "63% diff. max and eq. value", colour = "gray")+
  ggplot2::annotate("text", x = 75, y = 0.37*max(signal$y), label = "37% diff. max and eq. value", colour = "gray")+
  
  ggplot2::labs(x = "Time (arb. unit)",
           y = "Signal (arb. unit)",
           colour = "Legend")+
  ggplot2::theme(legend.position = "top", plot.title = ggplot2::element_text(hjust = 0.5))

## ----fig.width = 5, fig.height = 4, fig.align = "center", echo = FALSE--------
time <- 0:130
excitation <- c(rep(0,30),rep(1,50),rep(0,51))
signal <- generate.2order(time = time,
                            excitation = excitation,
                            y0 = 0,
                            v0 = 0,
                            xi = 0.2,
                            period = 15,
                            k = 1,
                            yeq = )
ggplot2::ggplot(signal)+
  ggplot2::ggtitle( "Second order differential equation solution")+
  ggplot2::geom_line(aes(t,y,color = "signal"))+
  ggplot2::geom_line(aes(time,excitation,color = "excitation"))+
  ggplot2::geom_vline(xintercept=80, colour = "gray")+
  ggplot2::annotate("text", x = 100, y = 1, label = "DLO", parse = TRUE, colour = "gray")+
  ggplot2::labs(x = "Time (arb. unit)",
           y = "Signal (arb. unit)",
           colour = "Legend")+
  ggplot2::theme(legend.position = "top", plot.title = ggplot2::element_text(hjust = 0.5))

