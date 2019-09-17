Canonical Correlation Analysis
==============================

Reference: http://www.nr.com/whp/notes/CanonCorrBySVD.pdf 

``` R code block (from cancor() help)
## signs of results are random
pop <- LifeCycleSavings[, 2:3]
oec <- LifeCycleSavings[, -(2:3)]
cancor(pop, oec)

x <- matrix(rnorm(150), 50, 3)
y <- matrix(rnorm(250), 50, 5)
(cxy <- cancor(x, y))
all(abs(cor(x %*% cxy$xcoef,
         y %*% cxy$ycoef)[,1:3] - diag(cxy $ cor)) < 1e-15)
all(abs(cor(x %*% cxy$xcoef) - diag(3)) < 1e-15)
all(abs(cor(y %*% cxy$ycoef) - diag(5)) < 1e-15)
```

``` Results of R cancor()
$cor
[1] 0.47379784 0.19096804 0.07742772

$xcoef
            [,1]        [,2]       [,3]
[1,] -0.03946325  0.10264493 0.12431097
[2,]  0.07220613 -0.08117339 0.08144780
[3,] -0.13680992 -0.07380634 0.03054464

$ycoef
             [,1]        [,2]         [,3]         [,4]         [,5]
[1,]  0.004633797 -0.15308796  0.011857994 -0.015677003 -0.009206857
[2,] -0.129339663  0.01571358  0.100287852  0.023048472  0.013950787
[3,]  0.041468863 -0.01611685  0.042866392  0.099784202 -0.010762694
[4,]  0.095474019  0.06161817  0.110423639 -0.069193087  0.010548204
[5,]  0.002062702 -0.04085351 -0.004747218 -0.008952294  0.120030530

$xcenter
[1] -0.12321636  0.04269298  0.11740652

$ycenter
[1] -0.20530476  0.19927215  0.08228493 -0.08209035  0.18698282

>      all(abs(cor(x %*% cxy$xcoef,
+                  y %*% cxy$ycoef)[,1:3] - diag(cxy $ cor)) < 1e-15)
[1] TRUE
>      all(abs(cor(x %*% cxy$xcoef) - diag(3)) < 1e-15)
[1] TRUE
>      all(abs(cor(y %*% cxy$ycoef) - diag(5)) < 1e-15)
[1] TRUE
```

``` Results of cancor.py
cor [0.8248 0.3653]
xcoef [[-0.0091 -0.0362]
 [ 0.0486 -0.2603]]
ycoef [[ 8.4710e-03  3.3379e-02 -5.1571e-03]
 [ 1.3074e-04 -7.5882e-05  4.5437e-06]
 [ 4.1706e-03 -1.2268e-02  5.1883e-02]]
xcenter pop15    35.0896
pop75     2.2930
dtype: float64
ycenter sr         9.6710
dpi     1106.7584
ddpi       3.7576
```
