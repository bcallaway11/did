# DRDID 1.0.3
  * Add new flags for non-unique unit identifier
  
  * Better handle of factor variables as covariates

# DRDID 1.0.2
  * Fix issue with NA in covariates
  
# DRDID 1.0.1
  * Allows for treating covariates as factor and alike when computing DiD
  
  * Improve error and warning handling due to collinearity and convergence issues.
  
# DRDID 1.0.0
  * First official version of package, functions for computing a variety of difference-in-differences (DiD) estimators for the ATT. 
  * Documentation is improved compared to the devel version, including examples for every function now.
  
  * Created wrapper function `drdid`, `ordid` and `ipwdid` to implement doubly-robust, outcome regression and inverse probability weighted DID estimators.
  
  * Add dataset used in the empirical application of Sant'Anna and Zhao (2020).



