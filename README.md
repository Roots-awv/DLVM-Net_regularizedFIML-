# A Regularized Maximum Likelihood Estimation of High-dimensional Network Models with Dynamic Latent Variables in the Context of Intensive Longitudinal Data
**Authors: Vivato V. Andriamiarana, Kevin Kistermann, Pascal Kilian & Augustin Kelava**
The manuscript is currently under revision in *Psychometrika*

## Repository Overview

This repository contains the R codes associated with the manuscript entitled Regularized Maximum Likelihood Estimation of High-Dimensional Network Models with Dynamic Latent Variables in the Context of Intensive Longitudinal Data.

The repository is organized into two main folders, corresponding to the empirical application and the simulation study.

### Empirical Application 
The first folder contains the R codes used to conduct the empirical analysis presented in the manuscript.

    - *Network_buidling_and_results.R* contains the code required to reproduce the empirical results reported in the manuscript. In particular, it covers the construction and visualization of the estimated network models, as well as the computation of confidence intervals for the parameter estimates.

    - *SAMdata12250-df17var7l1-FIML-aLassoFIML.R, SAMdata12250-df17var7l1-FIML-LassoFIML.R* and *SAMdata12250-df17var7l1-unregFIML.R* are illustrative implementations of the estimation approaches considered in the manuscript: the adaptive Lasso (aLasso), the Lasso, and the unregularized maximum likelihood approaches. Each script demonstrates the implementation of a specific estimation method using a particular penalty specification. Although the manuscript considers several penalty configurations and tuning parameter values, the repository includes only one illustrative example for each estimation approach.

    - *Computing_Penalized_Pseudo_FIML_LargeBT_CovMatrix.R* provides the computational functions required to implement the proposed estimation framework. In particular, it includes a collection of functions for matrix transformations and vector operations used to construct and manipulate the matrices involved in the estimation procedure. It also includes functions for covariance matrix operations, such as computing derivatives of the covariance matrix with respect to the model parameters, as well as functions for evaluating the full-information maximum likelihood (FIML) objective function and its gradient. In addition, this script contains functions for data management, descriptive statistics, and the generation of statistical summaries used in the empirical application.

### Simulation Study
The second folder contains the R scripts used to conduct the simulation study presented in the manuscript.

    - *mlDLVMNet-DGP.R* contains the code used to generate the simulated datasets according to the simulation design described in the manuscript.

    - *mlDLVMNet-adaptlasso-SimulationStudy.R, mlDLVMNet-lasso-SimulationStudy.R* and *mlDLVMNet-unreg-SimulationStudy.R* implement the estimation approaches considered in the simulation study: the adaptive Lasso (aLasso), the Lasso, and the unregularized maximum likelihood approaches. Each script contains the code required to estimate the model parameters using the corresponding estimation method.

    - *SimulationResults.R* provides the functions and procedures used to manage, organize, and process the simulation results.

