interference
============

MATLAB code for estimating crossover interference, following Campbell *et al*., 2014.

Requirements
------------

MATLAB 2014a or later, with Parallel Toolbox.

Usage
-----

The folder contains example functions to demonstrate simulation of crossover data given interference and escape paramters, followed by maximum likelihood estimation of those parameters.

+ **example_Stahl.m**      : Example of fitting the Housworth-Stahl model for simulated phase-known data.
+ **example_Stahl_quad.m** : Example of fitting the Housworth-Stahl model for simulated (phase-unknown) data.

The two scripts have similar usage syntax:

`[nu_est, p_est, lk_max, covariance_matrix] = example_Stahl(nu, p, cM_map_len, N_indv)`  
`[nu_est, p_est, lk_max, covariance_matrix] = example_Stahl_quad(nu, p, cM_map_len, N_indv)`

+ Output:
  - nu_est: estimated interference parameter.  
  - p_est: estimated escape parameter.  
  - lk_max: log likelihood of fitted model.  
  - covariance_matrix: covariance matrix estimated as inverse of fisher information matrix.  
  
+ Input:
  - nu: interference parameter to use for simulation (Required).  
  - p: escape parameter to use for simulation (Required).  
  - cM_map_len:	vector of chromosome map lengths in centiMorgans (optional. Default [200; 250]).  
  - N_indv: number of individuals to simulate (optional. Default 300).  
 
Other functions in the folder are used internally by the above scripts. Futher details can be found with the help comments at the start of each function. Brief descriptions follow:  

+ **simStahl.m**  : Simulate phase known data under the Housworth-Stahl model.
+ **simStahl_quad.m**  : Simulate phase unknown data (as derived from quartet families) under the Housworth-Stahl model.
+ **stahlLogLk.m** : Calculate the likelihood under the Housworth-Stahl model.
+ **stahlLogLk_quad.m** : Calculate the likelihood under the Housworth-Stahl model for phase-unknown data, as described in Campbell *et al*., 2014.
+ **derivest.m, fminsearchbnd.m, gradest.m, hessdiag.m, hessian.m** : All used internally, and derived directly from John D'Errico's contributions to the MATLAB file exchange (http://www.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd--fminsearchcon and http://www.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation).

