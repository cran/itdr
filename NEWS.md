---
output:
  pdf_document: default
  html_document: default
---
# *News*
# Changes in itdr 2.0.1 
* New Features

  - Integrate all the multivariate sufficient dimension reduction methods in regression within the _mitdr()_ function. 
  
  - Integrate _wh()_, _wx()_, and _wy()_ functions into a single function _hyperPara()_.


* Enhancements

  - Integrate all the uni-variate sufficient dimension reduction methods in regression within the _itdr()_ function. The _itdr()_ function now facilitate to use Fourier transformation method (FM), Convolution Transformation method (CM), Iterative hessian transformation method (iht), and inverse Fourier transformation method (invFM). 

  
* Bug Fixes

  - Fixed the errors in _recumbent_ dataset. 

# Changes in itdr 2.0.0 (2023-06-23)

* New Features

  - Include the following integral transformation methods.
  
    1). An iterated alternating direction method of multipliers (ADMM) algorithm that selects the sufficient variables using a Fourier transform sparse inverse regression estimators. This algorithm is integrated in _admmft()_  function. 
    
    2). A Minimum Discrepancy Approach with Fourier Transform in Sufficient Dimension Reduction. This algorithm is integrated in _fm_xire()_ function.  
  
* Enhancements

  - Include the following data sets to the package. 
 
    1). prostate - The data describe the level of a prostate-specific antigen associated with eight clinical measures in 97 male patients taking a radical prostatectomy.
    
    2). Raman - The Raman dataset contains 69 samples of fatty acid information in terms of percentage of total sample weight and percentage of total fat content
    
# Changes in itdr 1.2.0 (2022-04-12)

* Enhancements
  
  - Updated the package such that it matches with the R 4.2.0 upgrades. 


# itdr 1.0.0 (2021-07-23)

* **CRAN** Initial Submission.
* **GitHub** Initial Submission.







