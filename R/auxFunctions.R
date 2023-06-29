#-----------------------------------------------------------------------
# implements Eq. (1.5) of Munk and Tecuapetla-Gómez (2015)
dbe <- function(h, data) { 
  mean(diff(data, h)^2) / 2
}
#-----------------------------------------------------------------------
# implements Eq. (1.6) of Munk and Tecuapetla-Gómez (2015)
gammaZeroEst <- function(m, d, data) { 
  n <- length(data)
  sum((data[1:(n - 2 * (m + 1))] - (1 + d) * data[(m + 2):(n - m - 1)] + d * data[(2 * m + 3):n])^2) /
    (2 * (1 + d + d^2) * (n - 2 * (m + 1)))    
}
#
#----------------------------------------------------------------------- 
# implements Eq. (1.7) of Munk and Tecuapetla-Gómez (2015)
gh <- function(m, h, d, data) {
  ifelse(h <= m, gammaZeroEst(m=m, d=d, data=data) - dbe(h=h, data=data), 0)
} 
#

gh_first <- function(m, h, data) {
  ifelse(h <= m, dbe(h=m+1, data=data) - dbe(h=h, data=data), 0)
} 


#-----------------------------------------------------------------------
# computes first row of projection onto space of Toeplitz matrices
# see section 2.2.1 of Munk and Tecuapetla-Gómez (2015)
computeFirstRow <- function(matrix) {
  n <- ncol(matrix)
  firstRow <- numeric(n)

    for (i in 1:n) {
      for (j in 1:(n + 1 - i)) {
        firstRow[i] <- firstRow[i] + matrix[j, j - 1 + i]
      }
      firstRow[i] <- firstRow[i] / (n + 1 - i)
    }
  
  firstRow
}
#
#-----------------------------------------------------------------------
# tests a numeric for an integer value
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol  

# -----------------------------------------------------------------------------
# --- Added on June 23, 2023
dbacfSecondOrder <- function(data, m, d){
  
  if (missing(d)) {
    d <- numeric(m + 1L)
    d[1] <- 1
    if (m > 0) {
      for (i in 1:m) {
        aux <- i^2 - 4 * (m + 1 - i)^2
        if (aux >= 0) {
          d[i + 1L] <- (i + sqrt(aux)) / (2 * (m + 1 - i))
        } else {
          d[i + 1L] <- 1
        }
      }
    }      
  } else {
    if (length(d) != 1L && length(d) != m + 1L) {
      stop("length of d must be 1 or m + 1")
    }
    
    if ( !all(is.finite(d))) {
      stop("d must be a finite numeric")
    }
    
    if (length(d) == 1L) {
      d <- rep(d, m + 1L)
    }
  }
  
  acf <- numeric(m)
  if (m > 0) {
    for (h in 1:m)
      acf[h] <- gh(m=m, h=h, d=d[h + 1L], data=data)
  }
  acf <- c(gammaZeroEst(m=m, d=d[1], data=data), acf)   
  
  list(acf=acf, d=d)
}

# ----------------------------------------------------------------------------
# --- Added on June 23, 2023
dbacfFirstOrder <- function(data, m){
  acf <- numeric(m)
  if (m > 0) {
    for (h in 1:m)
      acf[h] <- gh_first(m=m+1, h=h, data=data) #dbe(k=h,data=data)
  }
  acf <- c(dbe(h=m+1, data=data), acf)    
  
  list(acf=acf, d=NULL)
}

