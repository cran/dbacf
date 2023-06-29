library(dbacf)

errors <- c(-0.1117,  0.8632,  0.4990,  0.8331,  0.1980,  2.0910,  0.3955,  0.4137,  0.1222,  1.2045)

n <- 10
m <- 1
gamma0 <- sum( (errors[1:(n - 2 *(m + 1))] - 2 * errors[(m + 2):(n - m - 1)] 
           + errors[(2 * (m + 1) + 1):n])^2 ) / (6 * (n - 2 * (m + 1)))

gamma1 <- gamma0 - sum((errors[1:(n - 1)] - errors[2:n])^2) / (2 * (n - 1))

acf.errors <- dbacf(errors, m = 1, order="second", plot = FALSE)$acf

stopifnot(all.equal(acf.errors, c(gamma0, gamma1), tol = 1e-12))