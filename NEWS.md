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
