# An R Implementation of Bayesian Optimization
Koji MAKIYAMA (@hoxo_m)  



[![Travis-CI Build Status](https://travis-ci.org/hoxo-m/bayesopt.svg?branch=master)](https://travis-ci.org/hoxo-m/bayesopt)
[![CRAN Version](http://www.r-pkg.org/badges/version/bayesopt)](http://cran.rstudio.com/web/packages/bayesopt)

## Usage

### Example 1. Sin Curve


```r
objective_func <- function(x) {
  x <- x * 10
  x * sin(x)
}

plot(objective_func)
```

![Figure 1: Sin Curve.](README_files/figure-html/unnamed-chunk-1-1.png)


```r
x <- seq(0, 1, length.out = 100)

library(bayesopt)
set.seed(314)
bo <- bayesopt(objective_func, x, iter = 10)
cat(sprintf("The best parameter is {'x': %f} with a score of %f", bo$opt_x, bo$opt_y))
```

```
## The best parameter is {'x': 0.818182} with a score of 7.746064
```

### Example 2. Himmelblau Function

![Figure 2: Himmelblau Function (from [Wikipedia](https://en.wikipedia.org/wiki/Himmelblau%27s_function)).](README_files/Himmelblau_Function.png)

The function has 4 answers.

- x = 3.0, y = 2.0, score = 0.
- x = -2.8, y = 3.1, score = 0.
- x = -3.8, y = -3.3, score = 0.
- x = 3.6, y = - 1.8, score = 0.


```r
x <- seq(0, 1, length.out = 100)
y <- seq(0, 1, length.out = 100)

trans_func <- function(x) (x - 0.5) * 12

Himmelblau <- function(x, y) {
  x <- trans_func(x)
  y <- trans_func(y)
  value <- (x^2 + y - 11)^2 + (x + y^2 - 7)^2
  -value
}

set.seed(314)
bo <- bayesopt(Himmelblau, x, y, iter = 20)
cat(sprintf("The best parameter is {'x': %f, 'y': %f} with a score of %f", trans_func(bo$opt_x[,1]), trans_func(bo$opt_x[,2]), bo$opt_y))
```

```
## The best parameter is {'x': 3.090909, 'y': 2.000000} with a score of -0.314869
```

## Related Work

- [rBayesianOptimization](https://cran.r-project.org/web/packages/rBayesianOptimization/index.html)

