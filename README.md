Canonical Correlation Analysis
==============================

Reference: http://www.nr.com/whp/notes/CanonCorrBySVD.pdf 

```R:R code block (from cancor() help)
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

```Results of R cancor()
$cor
[1] 0.8247966 0.3652762

$xcoef
              [,1]        [,2]
pop15 -0.009110856 -0.03622206
pop75  0.048647514 -0.26031158

$ycoef
             [,1]          [,2]          [,3]
sr   0.0084710221  3.337936e-02 -5.157130e-03
dpi  0.0001307398 -7.588232e-05  4.543705e-06
ddpi 0.0041706000 -1.226790e-02  5.188324e-02

$xcenter
  pop15   pop75
35.0896  2.2930

$ycenter
       sr       dpi      ddpi
   9.6710 1106.7584    3.7576

>
>      x <- matrix(rnorm(150), 50, 3)
>      y <- matrix(rnorm(250), 50, 5)
>      (cxy <- cancor(x, y))
$cor
[1] 0.3974386 0.2179173 0.1107121

$xcoef
            [,1]         [,2]        [,3]
[1,] -0.13609484 -0.005937912 -0.03212860
[2,] -0.02882604  0.126924567 -0.05283573
[3,]  0.06149630 -0.059445643 -0.10336198

$ycoef
             [,1]        [,2]        [,3]        [,4]        [,5]
[1,]  0.056676841  0.07036455  0.04587638 0.113444799  0.03066035
[2,] -0.115930714  0.12081643 -0.07117452 0.046568049 -0.02553839
[3,] -0.086182928 -0.01782606  0.07869779 0.007222671  0.01028472
[4,] -0.017922517 -0.09660660 -0.06424008 0.103308464 -0.02053949
[5,] -0.001728934 -0.01336417  0.03797456 0.002701664 -0.15885216

$xcenter
[1]  0.1630862 -0.1731617  0.2288064

$ycenter
[1]  0.03952011 -0.11775378 -0.01945354  0.24415264 -0.25356955

>      all(abs(cor(x %*% cxy$xcoef,
+                  y %*% cxy$ycoef)[,1:3] - diag(cxy $ cor)) < 1e-15)
[1] TRUE
>      all(abs(cor(x %*% cxy$xcoef) - diag(3)) < 1e-15)
[1] TRUE
>      all(abs(cor(y %*% cxy$ycoef) - diag(5)) < 1e-15)
[1] TRUE
```

```Results of cancor.py
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
