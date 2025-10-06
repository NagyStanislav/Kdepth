# SignDepth
R package *SignDepth*: Fast C++ and R implementation of the K-sign depth for univariate data, simplified 1-simplex depth and simplified 2-simplex depth for bivariate data, and the full 2-simplex depth for bivariate data.

To install the package *SignDepth* in Windows, you need to have *RTools* installed on your computer, see 

https://cran.r-project.org/bin/windows/Rtools/

Then it is enough to run

```R
#!R

library(devtools)
install_github("NagyStanislav/SignDepth")

library(SignDepth)
help(package="SignDepth")
```

All major functions in the package are implemented both in R and in C++ for computational efficacy.

- Christine H. Müller, Stanislav Nagy, and Samuel Trippler. (2025). Bivariate and multivariate sign depth and related distribution-free tests for model fit. _Under review._
