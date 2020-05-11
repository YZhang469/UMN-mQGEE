# Modified Q-learning with generalized estimating equations (GEE) for optimizing dynamic treatment regimes

This repository contains three R scripts for the Q-learning with GEE project and a text file of the data set used in the application section.

The R scripts display functions and syntax that are used to carry out the simulation study and data analysis. The results in the MQGEE paper can be reproduced by running the code in statistical software R. The following sections briefly explain how to use the files.

## Simulation study

### Custom functions

The document [functions.R](https://github.com/YZhang469/MQGEE/blob/master/functions.R) lists all the custom functions used in the simulation. It does not produce any results.

* The function `generateData.joint` generates a data set, where the repeated-measures outcomes follow a multivariate normal distribution. This is equivalent to a conditional data generating mechanism with a unmeasured covariate.
* The function `Q` analyzes the data set generated using composite Q-learning where the repeated-measures outcomes are collapsed with a weighted average. Users can choose between standard Q-learning and modified Q-learning. `Q` outputs a list of two models and two data sets: estimated models for stage 1 and stage 2 Q-functions using linear squares estimators, the original data set, and the estimated data set with optimal (counterfactual) outcomes based on parameter estimates from both Q-functions and the choice of standard or modified Q-learning.
* The function `QGEE` analyzes the data set generated using our proposed Q-learning with GEE, where the repeated-measures outcomes are treated as a vector. Likewise, users can choose between standard Q-learning and modified Q-learning. The function utilizes the R package [geepack](https://cran.r-project.org/web/packages/geepack/geepack.pdf) to perform estimation using GEE. Users have to specify the working correlations at both stage 1 and stage 2. The default is "unstructured" for stage 2 and "exchangeable" for stage 1. `QGEE` outputs a list of two models and two data sets: estimated stage 1 and stage 2 models using generalized estimating equations, the original data set, and the estimated data set with optimal (counterfactual) outcomes based on parameter estimated from both Q-functions and the choice of standard or modified Q-learning.

### Simulation code and results summary

The document [simulation.R](https://github.com/YZhang469/MQGEE/blob/master/simulation.R) displays the code for simulation studies performed in the modified Q-learning with GEE paper.

## Application

### Data set
To illustrate the implementation of modified Q-learning with GEE, we analyze a simulated data set (Nahum-Shani, 2020) from a simple version of ENGAGE study (McKay, 2015). The data set is readable from the text file [ENGAGEDataWideFormat.txt](https://github.com/YZhang469/MQGEE/blob/master/ENGAGEDataWideFormat.txt) and has 250 observations.

### Analysis

## References

McKay, J. R., Drapkin, M. L., Van Horn, D. H. A., Lynch, K. G., Oslin, D. W., De-Philippis, D., Ivey, M. and Cacciola, J. S. (2015). Effect of patient choice in an adaptivesequential randomization trial of treatment for alcohol and cocaine dependence. Journal of Consulting and Clinical Psychology, 83(6), 1021–1032.

Nahum-Shani, I., Almirall, D., Yap, J. R. T., McKay, J. R., Lynch, K. G., Freiheit,E. A. and Dziak, J. J.(2020). Smart longitudinal analysis: a tutorial for using repeatedoutcome measures from SMART studies to compare adaptive interventions. Psychological Methods, 25(1), 1–29.
