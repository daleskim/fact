# fact
Factor Analysis via Correlation Thresholding. An `R` package.

## Installation
Install the `devtools` package and run the following code:

`devtools::install_github("daleskim/fact")`

## Examples

```
# Factor Analysis Model
l <- .6
L <- matrix(c(rep(l, 3), rep(0, 4), rep(l, 3)), 5, 2)
P <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
O <- diag(diag(diag(5) - L %*% P %*% t(L)))
R <- L %*% P %*% t(L) + O
R
      [,1]  [,2]  [,3]  [,4]  [,5]
[1,] 1.000 0.360 0.468 0.108 0.108
[2,] 0.360 1.000 0.468 0.108 0.108
[3,] 0.468 0.468 1.000 0.468 0.468
[4,] 0.108 0.108 0.468 1.000 0.360
[5,] 0.108 0.108 0.468 0.360 1.000

# Run CT Algorithm (with specific thresholds)
ct_algorithm(R, .2)
[[1]]
     [,1] [,2]
[1,]    1    0
[2,]    1    0
[3,]    1    1
[4,]    0    1
[5,]    0    1

# Run CT Algorithm (with automatic thresholds)
ct_algorithm(R)
[[1]]
     [,1]
[1,]    1
[2,]    1
[3,]    1
[4,]    1
[5,]    1

[[2]]
     [,1] [,2]
[1,]    1    0
[2,]    1    0
[3,]    1    1
[4,]    0    1
[5,]    0    1

[[3]]
     [,1] [,2] [,3] [,4]
[1,]    1    0    0    0
[2,]    0    1    0    0
[3,]    1    1    1    1
[4,]    0    0    1    0
[5,]    0    0    0    1
```

## Reference
Kim, D. S., & Zhou, Q. (2023). Structure learning of latent factors via clique
search on correlation thresholded graphs. *Proceedings of the 40th
International Conference on Machine Learning*, 202, 16978–16996.
https://proceedings.mlr.press/v202/kim23aa.html
