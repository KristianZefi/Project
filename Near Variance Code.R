

library(expm)
library(Matrix)
library(matrixcalc)
library(mgcv)
Near_Var<-function (a,Convergence_Tolerance = 1e-10,
                    Iteration_Limit=300,Positive_eigen_Tolerance = 1e-10)
{W <- diag(1, nrow(a),ncol(a))
Projection_Onto_Set_S<-function(a){
  Wsqrt <- sqrtm(W)
  Winvsqrt <- matrix.inverse(Wsqrt)
  eigena<-eigen(a)
  eigenval<-eigena$values
  numbeig<-eigenval
  numbeig<-numbeig[numbeig>Positive_eigen_Tolerance]
  g<-length(numbeig)
  F<-slanczos(a,g,-1)
  brac<- F$vectors
  Eigen_Brac<-eigen(a)
  QBrac<-Eigen_Brac$vectors
  QValBrac<-Eigen_Brac$values
  D<-diag(QValBrac)
  D[D<Positive_eigen_Tolerance]<- Positive_eigen_Tolerance
  Spec_Brac<-QBrac %*% D %*% t(QBrac)
  PS<-Winvsqrt %*% Spec_Brac %*% Winvsqrt
  return(PS)
}
n <- nrow(a)
DS <- matrix(0,nrow(a),ncol(a))
Y<-a
j <- 0
converged <- FALSE
while (!converged ) {
  X <- Projection_Onto_Set_S(a)
  Y <- X
  Convergence_Test <- norm(Y - X, "I")/norm(Y, "I")
  j <- j + 1
  converged <- (Convergence_Test<= Convergence_Tolerance)
  if(j> Iteration_Limit){
    break }}
if (!converged)
  warning("DAPA did not converge")
eigen_Y <- eigen(Y, symmetric = TRUE)#Only using the
#positive Eiganvalues
eigenY_values <- eigen_Y$values
Positive_eigen_Tolerance <- Positive_eigen_Tolerance *abs(eigenY_values[1])
if (eigenY_values[n] < Positive_eigen_Tolerance) {
  eigenY_values[eigenY_values < Positive_eigen_Tolerance]<-Positive_eigen_Tolerance
  EigenY_Vectors <- eigen_Y$vectors
  Y <- EigenY_Vectors %*% (eigenY_values * t(EigenY_Vectors))
}
f<-eigen(Y)
EigenValues<-f$values
Rank<-qr(Y)$rank
structure(list(matrix=as.matrix(Y), Frobenius_Norm =
                 norm(a - Y, "F"),Rank = Rank,iterations =
                 j,eigenvalues=EigenValues))}

A<-matrix(c(5,3,2,3,5,5,2,5,5),nrow=3,ncol=3)
A
# 
# > A
# [,1] [,2] [,3]
# [1,]    5    3    2
# [2,]    3    5    5
# [3,]    2    5    5


# eigen(A)
# 
# > eigen(A)
# eigen() decomposition
# $values
# [1] 11.8390863  3.2893084 -0.1283947
# 
# $vectors
# [,1]       [,2]       [,3]
# [1,] -0.4612978  0.8691432  0.1783099
# [2,] -0.6465566 -0.1916741 -0.7383939
# [3,] -0.6075927 -0.4559069  0.6503692
# 

 Near_Var(A)
# 
# $matrix
# [,1]     [,2]     [,3]
# [1,] 5.004082 2.983095 2.014890
# [2,] 2.983095 5.070004 4.938341
# [3,] 2.014890 4.938341 5.054308
# 
# $Frobenius_Norm
# [1] 0.1283947
# 
# $Rank
# [1] 2
# 
# $iterations
# [1] 1
# 
# $eigenvalues
# [1] 1.183909e+01 3.289308e+00 1.183929e-09
# 
