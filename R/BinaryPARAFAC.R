#' Perform binary PARAFAC analysis on a three-dimensional dataset.
#'
#' This function implements binary PARAFAC analysis on a three-dimensional dataset,
#' using the user-specified algorithm.
#'
#' @param X The three-dimensional data matrix.
#' @param s The number of components to extract (default: 2).
#' @param I The number of levels in the first dimension of the data.
#' @param J The number of levels in the second dimension of the data.
#' @param K The number of levels in the third dimension of the data.
#' @param tolerance The tolerance for algorithm convergence (default: 1e-04).
#' @param penalization The value of the penalty to apply (default: 0.2).
#' @param num_max_iters The maximum number of allowed iterations (default: 100).
#' @param OptimMethod The optimization method to use (default: "CG").
#' @param Nombres The names of variables in the data matrix (optional).
#'
#' @return A list containing the results of binary PARAFAC analysis.
#'
#' @examples
#' BinaryPARAFAC(X = behaviour, s = 2, I = 128, J = 14, K = 11)
#'
#' @export
#'
#' @importFrom MultBiplotR RidgeBinaryLogistic
#' @importFrom stats optim pchisq
#' @importFrom rTensor khatri_rao
#'
BinaryPARAFAC <- function(X, s = 2, I, J , K, tolerance = 1e-04,  penalization=0.2,
                          num_max_iters=100, OptimMethod="CG",Nombres=NULL)
{
  b0=matrix(0,nrow=J*K, ncol=1)

  for (j in 1:(J*K))
    b0[j]=RidgeBinaryLogistic(y=X[,j], matrix(1,I,1), penalization = penalization)$beta

  A=matrix(1,nrow=I, ncol=1)
  B=matrix(1,nrow=J, ncol=1)
  C=matrix(1,nrow=K, ncol=1)
  CB= b0
  H=sigmoide(A %*% t(CB))
  L=sum(-1*X*log(H)-(1-X)*log(1-H), na.rm = TRUE)

  for (k in 1:s){
    parA=X[,k]
    parB=rep(1,J)/sqrt(J)
    parC=rep(1,K)/sqrt(K)
    A=cbind(A,parA)
    if (k == 1){
      B[,k]=parB
      C[,k]=parC
    }else
    {
      B = cbind(B,parB)
      C = cbind(C,parC)
    }
    err=1
    iter=0
    while( (err > tolerance) & (iter<num_max_iters)){
      iter=iter+1
      Lold=L

      #Update A
      resbipA <- optim(parA,fn=LLogRegARec,gr=grLogRegARec, method=OptimMethod, X=X, A=A, B= B, C= C, b = b0, lambda=penalization,r=k)
      parA=resbipA$par
      A[,k+1]=parA
      A[,k+1]=A[,k+1]-mean(A[,k+1])

      #Update B
      resbipB <- optim(parB, fn=LLogRegBRec, gr=grLogRegBRec, method=OptimMethod, X=X, A=A, B=B,C=C, b = b0,lambda=penalization,r=k)
      parB=resbipB$par
      B[,k]=parB

      #Update C
      resbipC <- optim(parC, fn=LLogRegCRec, gr=grLogRegCRec, method=OptimMethod, X=X, A=A, B=B,C=C,b = b0,lambda=penalization,r=k)
      parC=resbipC$par
      C[,k]=parC

      CB = cbind(b0,khatri_rao(C,B))
      H=sigmoide(A %*% t(CB))

      L=sum(-1*X*log(H)-(1-X)*log(1-H), na.rm = TRUE)
      err=abs(Lold-L)/L
      cat("\n",round(iter), round(L, 3), round(err,6))
    }
  }
  A=A[,-1]
  Res=list()
  Res$Data=X
  Res$Dimension=s
  Res$Penalization=penalization
  Res$Tolerance=tolerance
  Res$OptimMethod=OptimMethod
  Res$FirstMode=A
  Res$SecondMode=B
  Res$ThirdMode=C
  CB = cbind(b0,khatri_rao(C,B))
  esp = cbind(rep(1,I), A) %*% t(CB)
  pred = exp(esp)/(1 + exp(esp))
  pred2 = matrix(as.numeric(pred > 0.5), I, J*K)
  acier = matrix(as.numeric(round(X) == pred2), I, J*K)
  acierfil = 100*apply(acier,1,sum)/(J*K)
  aciercol = 100*apply(acier,2,sum)/I
  presences=apply(X, 2, sum)
  absences=I-presences
  sens = apply((acier==1) & (X==1), 2, sum)/presences
  spec = apply((acier==1) & (X==0), 2, sum)/absences
  totsens = sum((acier==1) & (X==1))/sum(presences)
  totspec = sum((acier==1) & (X==0))/sum(absences)
  gfit = (sum(sum(acier))/(I * J * K)) * 100


  Res$Biplot="Binary Logistic (PARAFAC)"
  Res$Type= "Binary Logistic (PARAFAC)"
  rownames(CB)=Nombres
  rownames(A) = rownames(X)
  Res$RowCoordinates=A
  Res$ColumnParameters=CB

  esp0 = matrix(rep(1,I), I,1) %*% CB[, 1]
  pred0 = exp(esp0)/(1 + exp(esp0))

  d1 = -2 * apply(X * log(pred0) + (1 - X) * log(1 - pred0),2,sum)
  d2 = -2 * apply(X * log(pred) + (1 - X) * log(1 - pred),2,sum)

  d = d1 - d2
  ps = matrix(0, J*K, 1)
  for (j in 1:(J*K)) ps[j] = 1 - pchisq(d[j], 1)

  Res$NullDeviances=d1
  Res$ModelDeviances=d2
  Res$Deviances=d
  Res$Dfs=rep(s, J*K)
  Res$pvalues=ps
  Res$CoxSnell=1-exp(-1*Res$Deviances/I)
  Res$Nagelkerke=Res$CoxSnell/(1-exp((Res$NullDeviances/(-2)))^(2/I))
  Res$MacFaden=1-(Res$ModelDeviances/Res$NullDeviances)

  Res$TotalPercent=gfit
  Res$ModelDevianceTotal=sum(Res$ModelDeviances)
  Res$NullDevianceTotal=sum(Res$NullDeviances)
  Res$DevianceTotal=sum(Res$Deviances)

  dd = sqrt(rowSums(cbind(1,Res$ColumnParameters[, 2:(s + 1)])^2))
  Res$Loadings = diag(1/dd) %*% Res$ColumnParameters[, 2:(s + 1)]
  Res$Tresholds = Res$ColumnParameters[, 1]/d
  Res$Communalities = rowSums(Res$Loadings^2)

  nn=I*J*K
  Res$TotCoxSnell=1-exp(-1*Res$DevianceTotal/nn)
  Res$TotNagelkerke=Res$TotCoxSnell/(1-exp((Res$NullDevianceTotal/(-2)))^(2/nn))
  Res$TotMacFaden=1-(Res$ModelDevianceTotal/Res$NullDevianceTotal)

  Res$R2 = apply((X-H)^2,2, sum)/apply((X)^2,2, sum)
  Res$TotR2 = sum((X-H)^2) /sum((X)^2)
  pred= matrix(as.numeric(H>0.5),I , J*K)
  verdad = matrix(as.numeric(X==pred),I , J*K)
  Res$PercentsCorrec=apply(verdad, 2, sum)/I
  Res$TotalPercent=sum(verdad)/(nn)
  Res$Sensitivity=sens
  Res$Specificity=spec
  Res$TotalSensitivity=totsens
  Res$TotalSpecificity=totspec
  Res$TotalDf = s*J*K
  Res$p=1-pchisq(Res$DevianceTotal, df = Res$TotalDf)

  Res$ClusterType="us"
  Res$Clusters = as.factor(matrix(1,I, 1))
  Res$ClusterColors="blue"
  Res$ClusterNames="ClusterTotal"
  class(Res) = "Binary.Logistic.Biplot"
}
