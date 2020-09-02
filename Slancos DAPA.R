library(expm)
library(matrix)
library(matrixcalc)
library(mgcv)

DAPA2<-function (a,Convergence_Tolerance = 1e-07,
                 Iteration_Limit=100,Positive_eigen_Tolerance = 1e-07)
{W <- diag(1, nrow(a),ncol(a))
Projection_Onto_Set_D<-function(a){
  diag(a)<-1
  return(a)}
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
  return(PS) }
n <- nrow(a)
DS <- matrix(0,nrow(a),ncol(a))
Y<-a
j <- 0
converged <- FALSE
while (!converged ) {
  R <- Y - DS
  X <- Projection_Onto_Set_S(R)
  DS<-X-R
  Y <- Projection_Onto_Set_D(X)
  Convergence_Test <- norm(Y - X, "I")/norm(Y, "I")
  j <- j + 1
  converged <- (Convergence_Test<= Convergence_Tolerance)
  if(j> Iteration_Limit){ break }}
if (!converged)
  warning("DAPA did not converge")
eigen_Y <- eigen(Y, symmetric = TRUE)
eigenY_values <- eigen_Y$values
Positive_eigen_Tolerance <- Positive_eigen_Tolerance *abs(eigenY_values[1])
if (eigenY_values[n] < Positive_eigen_Tolerance) {
  eigenY_values[eigenY_values < Positive_eigen_Tolerance]<- Positive_eigen_Tolerance
  EigenY_Vectors <- eigen_Y$vectors
  Y <- EigenY_Vectors %*% (eigenY_values * t(EigenY_Vectors))}
f<-eigen(Y)
EigenValues<-f$values
diag(Y) =1
Rank<-qr(Y)$rank
structure(list(as.matrix(Y), Frobenius_Norm =
                 norm(a - Y, "F"),Rank = Rank,iterations = j,eigenvalues=
                 EigenValues))}

A<-matrix(c(1,2,3,2,2,2,3,2,1),nrow=3,ncol=3)
A
eigen(A)
#eigen() decomposition
#$values
#[1]  6.000000e+00 -3.996803e-15 -2.000000e+00
DAPA2(A)
# #> DAPA2(A)
# #[[1]]
# "[,1]      [,2] [,3]
# [1,] 1.0000000 0.9999999    1
# [2,] 0.9999999 1.0000000    1
# [3,] 1.0000000 1.0000000    1
# 
# $Frobenius_Norm
# [1] 3.605551
# 
# $Rank
# [1] 1
# 
# $iterations
# [1] 39
# 
# $eigenvalues
# [1] 3e+00 3e-07 3e-07"

A<-matrix(c(1,1,0.2,1,1,0.295,0.2,0.295,1),nrow=3,ncol=3)
# A
# [,1]  [,2]  [,3]
# [1,]  1.0 1.000 0.200
# [2,]  1.0 1.000 0.295
# [3,]  0.2 0.295 1.000
# > 
eigen(A)
# > eigen(A)
# eigen() decomposition
# $values
# [1]  2.110531593  0.894250261 -0.004781854
# 
# $vectors
# [,1]       [,2]        [,3]
# [1,] -0.6674770  0.2633603  0.69650255
# [2,] -0.6810310  0.1623378 -0.71403307
# [3,] -0.3011166 -0.9509405  0.07099974

DAPA2(A)

# > DAPA2(A)
# [[1]]
# [,1]      [,2]      [,3]
# [1,] 1.0000000 0.9952874 0.2004666
# [2,] 0.9952874 1.0000000 0.2945218
# [3,] 0.2004666 0.2945218 1.0000000
# 
# $Frobenius_Norm
# [1] 0.006731326
# 
# $Rank
# [1] 3
# 
# $iterations
# [1] 15
# 
# $eigenvalues
# [1] 2.106240e+00 8.937599e-01 2.106240e-07



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
Near_Var(A)
