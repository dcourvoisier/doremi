### doremi v0.2.0
### New features
1. Added several methods for the estimation of derivatives in two-step procedures: functions calculate.glla, calculate.fda
2. calculate.gold handle non equidistant time points
2. Added input parameter "n" which is the maximum order of derivative to calculate (functions calculate.gold and calculate.glla)
3. Added function that find the optimum parameter for derivative estimation (embedding number in the case of gold and glla procedures and smoothing parameter in the case of fda).
Addition of a plot method to visualize the results of this function: evolution of R2 and parameter estimation with embedding number/smoothing parameter.
4. Added functions for the generation of simulation signals following a second order differential equation and the corresponding analysis tools for this type of signals
5. Modified the function that generates simulation signals. It now uses deSolve and it is able to generate a signal even for initial conditions different from t0=0
6. Summary from analysis function gives the coefficients estimted by the regression with their standard errors, together with the transformed coefficients characteristic from the model
7. Uses package futile.logger to log all events on the console (info, warning and errors). Info messages can be displayed by activating parameter "verbose=T"in analysis functions (first and second order).
8. Three vignettes are given: an introduction vignette, a first-order vignette, and the second order vignette

### doremi v0.1.1 
#### BUG FIXES
#### Vignette
1. Set seed in vignette so that the examples are reproductible.
#### Function remi
1. In single individual fit, time column in object "estimated" was fixed to "time" whereas now, it is the name given to the time variable in the input data frame, as done for the multiple subjects case.
2. When the data comes from a single individual and the parameter "id" is specified, there was an error because the function was trying to carry out a mixed-effects linear regression on a single individual.Corrected.
3. When specifying an incorrect id or signal parameter, the function didn't display the customized error coming from function errorcheck. It does now.
4. Function detects now if there are duplicated time points for each individual and displays an error message.
5. Bug when using multiple excitations and number of excitations was greater than 10, when renaming columns and numbering, numbers where passing from 9 to 20. This was due to the fact that variables are renamed to doremi fixed names and there was a problem with the implementation of gsub for obtaining original variable names again.

### doremi v0.1.0
First package release.
