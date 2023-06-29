library(dbacf)

# Toeplitz
b <- matrix(c(1, 1, 1, 1, 2, 1, 1, 1, 3), byrow = TRUE, nrow = 3)
bproj <- matrix(c(2, 1, 1, 1, 2, 1, 1, 1, 2), byrow = TRUE, nrow = 3)

c <- as.matrix(projectToeplitz(b))

stopifnot(identical(unname(bproj), unname(round(as.matrix(c), digits = 4))))

# positive semidefinite Toeplitz
d <- matrix(c(1, 1, 0, 1, 1, 1, 0, 1, 1), byrow = TRUE, ncol = 3)
dproj <- matrix(c(1, 0.7607, 0.1573, 0.7607, 1, 0.7607, 0.1573, 0.7607, 1), 
                byrow = TRUE, ncol = 3)

e <- nearPDToeplitz(d, type = "correlation")

stopifnot(identical(unname(dproj), unname(round(as.matrix(e$projection),
                                                digits = 4))))

