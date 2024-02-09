#' @importFrom rTensor khatri_rao
#' @importFrom ThreeWay permnew
#' @export
sigmoide<-function(z){
  (1/(1+exp(-1*z)))
}
#' @export
LLogRegARec <- function(par, X, A, B, C,b, lambda,r) { # Cost to estimate A
  A[,r+1]=par
  CB = cbind(b,khatri_rao(C,B))
  H=sigmoide(A %*% t(CB))
  L=sum(-1*X*log(H)-(1-X)*log(1-H), na.rm = TRUE)/2 + lambda*sum(A^2, na.rm = TRUE)/2
  return(L)
}
#' @export
LLogRegBRec <- function(par, X, A, B, C, b, lambda,r) { # Cost to estimate B
  B[,r]=par

  I = dim(A)[1]
  J = dim(B)[1]
  K = dim(C)[1]
  CB = cbind(b,khatri_rao(C,B))
  H=sigmoide(A %*% t(CB))
  Xb=permnew(X,I,J,K)
  Hb=permnew(H,I,J,K)
  L=sum(-1*Xb*log(Hb)-(1-Xb)*log(1-Hb), na.rm = TRUE)/2 + lambda*sum(B^2, na.rm = TRUE)/2
  return(L)
}
#' @export
LLogRegCRec <- function(par, X, A, B, C, b, lambda,r) { # Cost to estimate C
  C[,r] = par
  I = dim(A)[1]
  J = dim(B)[1]
  K = dim(C)[1]
  CB = cbind(b,khatri_rao(C,B))
  H=sigmoide(A %*% t(CB))
  Xc=permnew(X,J,K,I)
  Hc=permnew(H,J,K,I)
  L=sum(-1*Xc*log(Hc)-(1-Xc)*log(1-Hc), na.rm = TRUE)/2 + lambda*sum(C^2, na.rm = TRUE)/2
  return(L)
}

#' @export
grLogRegARec <- function(par, X, A, B, C, b, lambda,r) { ## Gradient to estimate A
  A[,r+1]=par
  CB = cbind(b,khatri_rao(C,B))
  H=sigmoide(A %*% t(CB))
  E = H-X
  E[which(is.na(X))]=0
  gradA=E%*%khatri_rao(C,B)+lambda*A[,-1]
  grad=c(c(gradA[,r]))
  return(grad)
}
#' @export
grLogRegBRec <- function(par, X,A, B, C,b, lambda,r) { ## Gradient to estimate B
  B[,r]=par
  I = dim(A)[1]
  J = dim(B)[1]
  K = dim(C)[1]
  CB = cbind(b,khatri_rao(C,B))
  H=sigmoide(A %*% t(CB))
  Xb=permnew(X,I,J,K)
  Hb=permnew(H,I,J,K)
  E = Hb-Xb
  E[which(is.na(Xb))]=0
  AC = khatri_rao(as.matrix(A[,-1]),C)
  gradB=E%*%AC+lambda*B
  grad=c(c(gradB[,r]))
  return(grad)
}
#' @export
grLogRegCRec <- function(par, X,A, B, C, b, lambda,r) { ## Gradient to estimate C
  C[,r]=par
  I = dim(A)[1]
  J = dim(B)[1]
  K = dim(C)[1]
  CB = cbind(b,khatri_rao(C,B))
  H=sigmoide(A %*% t(CB))
  Xc=permnew(X,J,K,I)
  Hc=permnew(H,J,K,I)
  E = Hc-Xc
  E[which(is.na(Xc))]=0
  BA = khatri_rao(B,as.matrix(A[,-1]))
  gradC=E%*%BA+lambda*C
  grad=c(c(gradC[,r]))
  return(grad)
}
