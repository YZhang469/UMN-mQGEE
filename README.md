# Modified Q-learning with generalized estimating equations (GEE) for optimizing dynamic treatment regimes

This repository contains three R scripts for the Q-learning with GEE project. The R scripts display functions and syntax that are used to carry out the simulation study and data analysis. The results in the MQGEE paper can be reproduced by running the code in statistical software R. The following sections briefly explain how to use the files.

## Simulation study

### Custom functions

The document [functions.R](https://github.com/YZhang469/MQGEE/blob/master/functions.R) lists all the custom functions used in the simulation. It does not produce any results.

* The function `generateData.joint` generates a data set, where the repeated-measures outcomes follow a multivariate normal distribution. This is equivalent to a conditional data generating mechanism with a unmeasured covariate.
* The function `Q` analyzes the data set generated using composite Q-learning where the repeated-measures outcomes are collapsed with a weighted average. Users can choose between composite Q-learning (SQ) and modified composite Q-learning (MQ). `Q` outputs a list of two models and two data sets: estimated models for stage 1 and stage 2 Q-functions using linear squares estimators, the original data set, and the estimated data set with optimal (counterfactual) outcomes based on parameter estimates from both Q-functions and the choice of standard or modified composite Q-learning.
* The function `QGEE` analyzes the data set generated using our proposed Q-learning with GEE, where the repeated-measures outcomes are treated as a vector. Likewise, users can choose between GEE Q-learning (SQGEE) and modified GEE Q-learning (MQGEE). The function utilizes the R package [geepack](https://cran.r-project.org/web/packages/geepack/geepack.pdf) to perform estimation using GEE. Users have to specify the working correlations at both stage 1 and stage 2; the default is "unstructured" for stage 2 and "exchangeable" for stage 1. `QGEE` outputs a list of two models and two data sets: estimated stage 1 and stage 2 models using generalized estimating equations, the original data set, and the estimated data set with optimal (counterfactual) outcomes based on parameter estimated from both Q-functions and the choice of standard or modified GEE Q-learning.

### Simulation code and results

The document [simulation.R](https://github.com/YZhang469/MQGEE/blob/master/simulation.R) contains the code for simulation studies performed in the MQGEE paper.

* The function `metric` defines the metrics we use to compare the performance of different methods. It returns a list of four elements for an underlying data set: `pci.s2`, a scalar value indicating the empirical probability of correctly identifying the stage 2 optimal rule; `pci.s1`, a scalar value indicating the empirical probablity of correctly identifying the stage 1 optimal rule; `rmse`, a vector of three elements, with each element being the root mean square error of the estimates of heterogeneous (individual/conditional) causal effect of at each time point; `bias.mat`, a 3\*7 matrix, where each column represents the bias of causal effect estimates at various time points based on each value of the baseline measurement Y<sub>0</sub> in {-3, -2, -1, 0, 1, 2, 3}.

* The function `simulation` runs I simulations, where I is the number of simulations and is user-specified, and generates a list of metrics associated with the simulated data sets.

* The rest of the script shows the code to perform the simulation under different scenarios discussed in the MQGEE paper, and summarize the results using tables and figures. Note that it usually takes from 30 minutes to 1 hour for 1000 simulations to run.

## Application

### Data set

To illustrate the implementation of modified GEE Q-learning, we analyze a simulated data set (Nahum-Shani, 2020) from a simple version of ENGAGE study (McKay, 2015). ENGAGE is a sMART study for patients suffering from relatively severe substance-use disorders. The data set is readable from [5.TutorialDatasetWideFormat.txt](http://supp.apa.org/psycarticles/supplemental/met0000219/5.TutorialDatasetWideFormat.txt).

The format of the data is (X<sub>0</sub>, A<sub>1</sub>, Y<sub>1</sub>, R, A<sub>2</sub>, Y<sub>2</sub>, Y<sub>3</sub>) with a sample size of 250. X<sub>0</sub> represents the baseline covariates including age, gender, an indicator of post-secondary education, and the number of heavy drinking days. A<sub>1</sub> is a binary variable with values from {-1, 1}, indicating the treatment received at stage 1. Y<sub>1</sub> is the treatment readiness, which is the outcome of interest, measured right after stage 1 treatment period. R is a tailoring variable indicating whether an individual responded to the assigned stage 1 treatment. A<sub>2</sub> takes values from {-1, 0, 1}, depending on the embedding tailoring. Y<sub>2</sub> and Y<sub>3</sub> are the repeated-measures outcome at the end of stage 2.

### Analysis

The function `analyzeData` analyzes the abovementioned data set using Q-learning with GEE in the context of the simplified ENGAGE study. Users can specify whether they want to use standard Q-learning or modified Q-learning. Results are then summarized and relevant figures are plotted.

## References

McKay, J. R., Drapkin, M. L., Van Horn, D. H. A., Lynch, K. G., Oslin, D. W., De-Philippis, D., Ivey, M. and Cacciola, J. S. (2015). Effect of patient choice in an adaptivesequential randomization trial of treatment for alcohol and cocaine dependence. Journal of Consulting and Clinical Psychology, 83(6), 1021–1032.

Nahum-Shani, I., Almirall, D., Yap, J. R. T., McKay, J. R., Lynch, K. G., Freiheit,E. A. and Dziak, J. J.(2020). Smart longitudinal analysis: a tutorial for using repeatedoutcome measures from SMART studies to compare adaptive interventions. Psychological Methods, 25(1), 1–29.
