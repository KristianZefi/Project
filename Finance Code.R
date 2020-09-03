MinVarPortfolio<-function(estvar){
  s<-diag(estvar)
  n<-nrow(estvar)
  c<-ncol(estvar)
  if(c!=n){stop("not a square Covariance Matrix")}
  if(any(s<0))
  {stop("Covariance matrix seems to be negative semi-definite")}
  else {
    vect1<-as.matrix(rep(1,each=n))
    Numerator<-solve(as.matrix(estvar),tol=1e-20)%*%vect1
    denom1<-t(vect1) %*%matrix.inverse(as.matrix(estvar))
    denom2<-denom1 %*%vect1
    f<-Numerator/denom2[1,1]}
  print(f)
}

# 
# > B<-matrix(c(135.68,131.74,137.95,138.05 ,142.53,144.23,1539.13 ,1500.28,1575.39 ,1629.51 ,1656.58,1659.42,136.22 
#,127.99  ,136.19 ,143.40 ,139.83,142.58,101.12 ,97.40,101.93,102.06 ,102.80 ,104.27,109.28,108.26,110.48,109.53,110.58
#,107.49),6,5)
# > B
#stock prices of Facebook, Amazon, Nvidia, Microsoft, Pepsi
# [,1]    [,2]   [,3]   [,4]   [,5]
# [1,] 135.68 1539.13 136.22 101.12 109.28
# [2,] 131.74 1500.28 127.99  97.40 108.26
# [3,] 137.95 1575.39 136.19 101.93 110.48
# [4,] 138.05 1629.51 143.40 102.06 109.53
# [5,] 142.53 1656.58 139.83 102.80 110.58
# [6,] 144.23 1659.42 142.58 104.27 107.49

# 
# 
# > c<-cov(B)
# > MinVarPortfolio(c)
# Global Minimum Variance Portfolio of the 5 stocks
# [1,]  1.40112788
# [2,] -0.09072429
# [3,]  0.66989498
# [4,] -1.71849855
# [5,]  0.73819999


# 
# > A
#      [,1] [,2] [,3]
# [1,]    5    3    2
# [2,]    3    5    5
# [3,]    2    5    5
# > b<-cov(A)
# > MinVarPortfolio(b)
# Error in solve.default(x) : 
#   system is computationally singular:
#   reciprocal condition number = 1.77636e-17 
##Matrix Inverse is Computationally impossible to solve so implement Pseudo Inverse



PseudInv<-function(B){
  eig<-eigen(B)
  val<-eig$values
  svdb<-svd(B)
  U<-svdb$u
  D<-diag(1/svdb$d)
  D[,which(val<1e-05)]<-1e-05
  v<-svdb$v
  a<-v%*%D
  a%*%t(U)
}

MinVarPortfolio_Pseudo<-function(estvar){
  s<-diag(estvar)
  n<-nrow(estvar)
  c<-ncol(estvar)
  if(c!=n){stop("not a square Covariance Matrix")}
  if(any(s<0))
  {stop("Covariance matrix seems to be negative semi-definite")}
  else {
    vect1<-as.matrix(rep(1,each=n))
    Numerator<-solve(as.matrix(estvar),tol=1e-20)%*%vect1
    denom1<-t(vect1) %*%PseudInv(as.matrix(estvar))
    denom2<-denom1 %*%vect1
    f<-Numerator/denom2[1,1]}
  print(f)
}

# > MinVarPortfolio_Pseudo(b)
# [,1]
# [1,]  4.294469e-01
# [2,] -9.480671e+13
# [3,]  6.320447e+13
