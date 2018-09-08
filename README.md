LEP
=======

LEP is a statistical approach to integrating individual level genotype data and summary statistics by LEveraging Pleiotropy. 'LEP' R package provides computationally efficient and user friendly interface to fit and evaluate the LEP model. It accepts both the R-type data  and binary plink files.

Usage
=======

[The 'LEP' vignette](https://github.com/daviddaigithub/LEP/blob/master/inst/doc/LEP_package.pdf) will provide a good start point for the genetic analysis using LEP package.


Development 
=======
This R package is developed by Mingwei Dai and Can Yang, and maintained by Can Yang <eeyangc@gmail.com>.

Installation
=======
To install the development version of LEP, it's easiest to use the 'devtools' package. Note that LEP depends on the 'Rcpp' and 'RcppArmadillo' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.

install.packages("devtools")  
library(devtools)  
install_github("daviddaigithub/LEP")  

References
=======
M. Dai, X. Wan, H. Peng, Y. Wang, Y. Liu, J. Liu, Z. Xu and C. Yang. LEP: Joint Analysis of Individual Level Genotype Data and Summary Statistics by Leveraging Pleiotropy. Submitted


