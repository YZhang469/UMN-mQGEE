# R Code and example data for "Modified Q-learning with generalized estimating equations for optimizing dynamic treatment regimes with repeated-measures outcomes"

This file contains five R scripts which display functions and syntax that are used to carry out the simulation study and data analysis for the modified Q-learning with GEE (mQGEE) paper. The results in the mQGEE paper can be reproduced by running the code in statistical software R version 4.0.3 (2020-10-10). The following sections briefly explain how to use the files.

## Simulation study

### Custom functions

The document [functions.R](https://github.com/YZhang469/MQGEE/blob/master/functions.R) lists all the custom functions used in the simulation.

* The function `generateData.joint` generates a dataset, where the repeated-measures outcomes follow a multivariate normal distribution. This is equivalent to a conditional data generating mechanism with an unmeasured covariate.
* The function `Q` analyzes the dataset generated using composite Q-learning where the repeated-measures outcomes are collapsed with a weighted average. Users can choose between composite Q-learning (Q) and modified composite Q-learning (mQ). `Q` outputs a list of two models and two datasets: estimated models for stage 1 and stage 2 Q-functions using linear squares estimators, the original dataset, and the estimated dataset with optimal (counterfactual) outcomes based on parameter estimates from both Q-functions.
* The function `QGEE` analyzes the dataset generated using our proposed Q-learning with GEE, where the repeated-measures outcomes are treated as a vector. Likewise, users can choose between GEE Q-learning (QGEE) and modified GEE Q-learning (mQGEE). The function utilizes the R package [geepack](https://cran.r-project.org/web/packages/geepack/geepack.pdf) to perform estimation using GEE. Users have to specify the working correlations at both stage 1 and stage 2; the default is "unstructured" for stage 2 and "exchangeable" for stage 1. `QGEE` outputs a list of two models and two datasets: estimated stage 1 and stage 2 models using generalized estimating equations, the original dataset, and the estimated dataset with optimal (counterfactual) outcomes based on parameter estimated from both Q-functions.
* The function `metric` defines the metrics we use to compare the performance of different methods. It returns a list of four elements for an underlying dataset: `pci.s2`, a scalar value indicating the empirical probability of correctly identifying the stage 2 optimal rule; `pci.s1`, a scalar value indicating the empirical probablity of correctly identifying the stage 1 optimal rule; `rmse`, a vector of three elements, with each element being the root mean square error of the estimates of heterogeneous (individual/conditional) causal effect of at each time point; `bias.mat`, a 3\*7 matrix, where each column represents the bias of causal effect estimates at various time points based on each value of the baseline measurement Y<sub>0</sub> in {-3, -2, -1, 0, 1, 2, 3}.

### Simulation code and results

The document [simulation_by_sample_size.R](https://github.com/YZhang469/UMN-mQGEE/blob/master/simulation_by_sample_size.R) shows the code for the simulation study in the mQGEE manuscript.

* The function `simulation` runs *I* simulations, where *I* is the number of simulations and is user-specified, and generates a list of metrics associated with the simulated datasets.

* The rest of the script shows the code to perform the simulation under different scenarios discussed in the mQGEE paper, and summarize the results using tables and figures. Note that it usually takes from 30 minutes to 1 hour for 1000 simulations to run.

The documents [simulation_by_weights.R](https://github.com/YZhang469/UMN-mQGEE/blob/master/simulation_by_weights.R) and [simulation_relative_efficiency.R](https://github.com/YZhang469/UMN-mQGEE/blob/master/simulation_relative_efficiency.R) contain the code for the simulation studies in the mQGEE supplementary materials.

## Application

### Data

To illustrate the implementation of mQGEE, we analyze a dataset collected from the M-bridge study. While we are not able to provide the original dataset, we include a simulated dataset [simulated_data.csv](https://github.com/YZhang469/UMN-mQGEE/blob/master/simulated_data.csv) which has the same format as the M-bridge data.

### Analysis

The document [data_analysis.R](https://github.com/YZhang469/UMN-mQGEE/blob/master/data_analysis.R) shows the code to analyze the abovementioned dataset

* The function `analyzeData` analyzes the dataset using mQGEE. Users can specify whether they want to use QGEE or mQGEE. Since we have a considerably large set of covariates (with time-dependent coefficients) in the model at each stage, we perform variable selection using the function `gee_stepper`. Results are then summarized using relevant tables and figures.

* The R code for the additional results in the supplementary materials, where unequal weights of the longitudinal outcomes are considered, is also included in the document.
