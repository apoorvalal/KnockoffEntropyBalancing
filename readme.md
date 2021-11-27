Entropy Balancing in High Dimensions
====================================

**Apoorva Lal**

Minimal working example of Double-LASSO augmented entropy balancing
estimator proposed in Lal (2021).

``` r
rm(list = ls())
library(LalRUtils)
set.seed(42)
# paths
source("estimators.R")
```

    ##       wants        loaded
    ##  [1,] "tidyverse"  TRUE  
    ##  [2,] "data.table" TRUE  
    ##  [3,] "hdm"        TRUE  
    ##  [4,] "ebal"       TRUE  
    ##  [5,] "glue"       TRUE  
    ##  [6,] "rlist"      TRUE  
    ##  [7,] "lfe"        TRUE  
    ##  [8,] "janitor"    TRUE  
    ##  [9,] "caret"      TRUE  
    ## [10,] "kosel"      TRUE

experimental benchmark (Lalonde 86)
---------------------------

``` r
library(causalsens)
data(lalonde.exp)
summary(robustify(felm(re78 ~ treat, lalonde.exp)))$coefficients[2, 1:2]
```

    ##   Estimate Std. Error 
    ##     1794.3      670.8

Implementation on PSID Lalonde Sample
-------------------------------------

Estimation akin to Dehejia & Wahba (1999).

### Traditional Estimators

``` r
data(lalonde.psid)
dtpsid = data.table(lalonde.psid)
y = 're78'
w =  'treat'
x = setdiff(colnames(lalonde.psid), c(y, w))
# %% traditional estimators
dm_est  = diffmeans(dtpsid, y, w)
cov_est = covaradjust(dtpsid, y, w, x)
ipw_est = ipw_reg(dtpsid, y, w, x)
```

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

``` r
eb_est = eb_reg(dtpsid, y, w, x)
lr_est = rbind(dm_est, cov_est, ipw_est, eb_est)
rownames(lr_est) = c('difference in means', 'regression', 'ipw', 'entropy balancing')
```

### Double LASSO and Knockoff EB

``` r
y = 're78'
w =  'treat'
x = setdiff(colnames(lalonde.psid), c(y, w))
# %% basis expansion and interaction creation step
data = prep_matrices('re78', 'treat', controls = x, dat = dtpsid)
y = data[[1]]
w = data[[2]]
X = data[[3]]

# %% estimation using double lasso, and knockoff lasso + knockoff-ebal
hdres = estimate_hdc(y, w, X)
```

final estimates
---------------

``` r
rbind(lr_est, hdres)
```

    ##                     Estimate Std. Error
    ## difference in means   -15205      655.9
    ## regression             -8746     3720.6
    ## ipw                   -14501     1132.9
    ## entropy balancing       2425      876.5
    ## Partial-Out            -1553     1141.1
    ## Double Selection       -1651      938.3
    ## Knockoff Selector      -3717      916.7
    ## Knockoff EB             2302      945.3

Knockoff-EB comes closest to the experimental benchmark. All the
regression estimators (and IPW) estimate the wrong sign.
