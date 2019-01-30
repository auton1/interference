interference
============

Python implementation of code for estimating crossover interference under the Housworth-Stahl model, following Campbell *et al*., 2014.

Requirements
------------

Python 2.7, numpy, acipy, etc.

Usage
-----

The folder contains example functions to demonstrate simulation of crossover data given interference and escape paramters, followed by maximum likelihood estimation of those parameters.

+ **example_Stahl.py**      : Example of fitting the Housworth-Stahl model for simulated data (both phase-known and phase-unknown).

The two scripts have similar usage syntax:

`[nu_est, p_est] = example_Stahl(nu, p, cM_map_len, n_indv)`  
`[nu_est, p_est] = example_Stahl_quad(nu, p, cM_map_len, n_indv)`

+ Output:
  - nu_est: estimated interference parameter.  
  - p_est: estimated escape parameter.  
  
+ Input:
  - nu: interference parameter to use for simulation (Required).  
  - p: escape parameter to use for simulation (Required).  
  - cM_map_len:	vector of chromosome map lengths in centiMorgans (optional. Default [200; 250]).  
  - n_indv: number of individuals to simulate (optional. Default 300).  
 
Other functions in the folder are used internally by the above scripts. Futher details can be found with the help comments at the start of each function. Brief descriptions follow:  

+ **simStahl.py**  : Simulate data under the Housworth-Stahl model (both phase-known and phase-unknown routines).
+ **stahlLogLk.py** : Calculate the likelihood under the Housworth-Stahl model for phase-known data, and phase-unknown data, as described in Campbell *et al*., 2014.

Performance
-----------

The code has been tested against the corresponding MATLAB implementation. 
