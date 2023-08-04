# R1
l <- .6
L <- matrix(c(rep(l, 6), rep(0, 6), rep(l, 6)), 9, 2)
P <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
O <- diag(diag(diag(9) - L %*% P %*% t(L)))
R1 <- L %*% P %*% t(L) + O
tau <- seq(0, 1, length.out = 50)
rm(l, L, P, O)

# R2 - Randomly Negated R1
R2 <- matrix(0, 9, 9)
R2[lower.tri(R2)] <- R1[lower.tri(R1)] * sample(c(-1, 1), 36, replace = TRUE)
R2[upper.tri(R2)] <- t(R2)[upper.tri(R2)]


