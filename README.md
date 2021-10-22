# NetGAM
### By: Samantha J. Gleich, Jacob A. Cram, Jake L. Weissman, and David A. Caron
### Last Modified: October 17, 2021
`NetGAM` is an R package that is designed to decrease the influence that time-series properties (i.e. seasonality, long-term trends, and autocorrelation) have on ecological network analyses. The `NetGAM` package takes species abundance data as input and uses generalized additive mixed models (gamms) in the `mgcv` package to model each species as a function of time-series predictor variables. The residuals of each gamm are then extracted and used in downstream networking methods (e.g. graphical lasso, correlation networks, etc.).

The `NetGAM` package requires the use of other R packages:
- `mgcv`
- `dplyr`
- `compositions`
- `stats`
- `psych`
- `pulsar`
- `batchtools`
- `huge`

# Installation
The `NetGAM` package can be installed using the `devtools` package
```
install.packages("devtools")
library(devtools)
devtools::install_github("sgleich/NetGAM")
library(NetGAM)
```
The `NetGAM` package dependencies can also be installed as follows:
```
install.packages(c("mgcv","dplyr","compositions","stats","psych","pulsar","batchtools","huge"))
```

# Example
An example of how to use the NetGAM package is outlined here file:///private/var/folders/5n/08h315tn3n12xmlntsd5w2s80000gn/T/RtmpiroxBh/preview-17d7433bb14a2.dir/Example_NetGAM_Network.html
