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