# A Regularized Maximum Likelihood Estimation of High-dimensional Network Models with Dynamic Latent Variables in the Context of Intensive Longitudinal Data
**Authors: Vivato V. Andriamiarana, Kevin Kistermann, Pascal Kilian & Augustin Kelava**

The manuscript is currently under revision in *Psychometrika*.

## Repository Overview

This repository contains the R codes associated with the manuscript entitled Regularized Maximum Likelihood Estimation of High-Dimensional Network Models with Dynamic Latent Variables in the Context of Intensive Longitudinal Data.

The repository is organized into two main folders, corresponding to the empirical application and the simulation study.

### Main functions
- `Computing_Penalized_Pseudo_FIML_LargeBT_CovMatrix.R` provides the computational functions required to implement the proposed estimation framework. In particular, it includes a collection of functions which can be categiruzed as follow:
    - *Transformation matrix/vector* which contains matrix and vector operations used to construct and manipulate the matrices involved in the estimation procedure.  
    - *Covariance matrices structure* which includes functions for covariance matrix structures
    - *The full information likelihood function* which is the objective function required for `optim()` and its "L-BFGS-B" algorithm.
    - *The score function* which includes the gradinet function required for `optim()` and olso the derivatives of the covariance matrix with respect to the model parameters.
    - *Data management* which contains a function used for data management for time series, the sample means and the sample-implied covariance matrix.

### Empirical Application 
The first folder contains the R codes used to conduct the empirical analysis presented in the manuscript.

- `Network_buidling_and_results.R` contains the code used for handling the empirical results reported in the manuscript. In particular, it covers the construction and visualization of the estimated network models, as well as the computation of confidence intervals for the parameter estimates.

- The following syntaxes correspind estimation methods considered in the manuscript: 
    - `SAMdata12250-df17var7l1-FIML-aLassoFIML.R` correponds to the syntax of the adaptive Lasso (aLasso); 
    - `SAMdata12250-df17var7l1-FIML-LassoFIML.R` corresponds to the syntax for the Lasso;
    - `SAMdata12250-df17var7l1-unregFIML.R` corresponds to the unregularized approach. 

### Simulation Study
The second folder contains the R scripts used to conduct the simulation study presented in the manuscript.

- `mlDLVMNet-DGP.R` contains the code used to generate the simulated datasets according to the simulation design described in the manuscript.

- The simulation study consider three estimaton approaches for which the respective syntaxes are: 
    - `mlDLVMNet-adaptlasso-SimulationStudy.R` implementing the adaptive Lasso approaches
    - `mlDLVMNet-lasso-SimulationStudy.R` implementing the Lasso approaches
    - `mlDLVMNet-unreg-SimulationStudy.R` which implements the unregularized approach.

- `SimulationResults.R` provides the functions and procedures used to manage, organize, and process the simulation results.

