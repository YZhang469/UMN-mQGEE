# Modified Q-learning with generalized estimating equations (GEE) for optimizing dynamic treatment regimes

This repository contains three R scripts for the Q-learning with GEE project and a text file of the data set used in the application section.

The R scripts display functions and syntax that are used to carry out the simulation study and data analysis. The results in the MQGEE paper can be reproduced by running the code in R. The following sections briefly explain how to use the files.

## Custom functions

The document functions.R lists all custom functions for Q-learning used in the simulation. It does not produce any results.

* The function `generateData.joint`

## Simulation code

The document simulation.R displays the code for simulation studies performed in the modified Q-learning with GEE paper.

## Application

To illustrate the implementation of modified Q-learning with GEE, we analyze a simulated data set (Nahum-Shani, 2020) from a simple version of ENGAGE study (McKay, 2015). The data set is readable from the text file ENGAGEDataWideFormat.txt and has 250 observations.

### References

McKay, J. R., Drapkin, M. L., Van Horn, D. H. A., Lynch, K. G., Oslin, D. W., De-Philippis, D., Ivey, M. and Cacciola, J. S. (2015). Effect of patient choice in an adaptivesequential randomization trial of treatment for alcohol and cocaine dependence. Journal of Consulting and Clinical Psychology, 83(6), 1021–1032.

Nahum-Shani, I., Almirall, D., Yap, J. R. T., McKay, J. R., Lynch, K. G., Freiheit,E. A. and Dziak, J. J.(2020). Smart longitudinal analysis: a tutorial for using repeatedoutcome measures from SMART studies to compare adaptive interventions. Psychological Methods, 25(1), 1–29.
