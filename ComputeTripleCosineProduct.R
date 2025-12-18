ComputeTripleCosineProduct <- function( A , B , C )
{
  #A = Sol.A
  #B = Sol.B
  #C = Sol.C
  nComp = dim(A)[2]
  outA = matrix( NA , nComp , nComp )
  outB = matrix( NA , nComp , nComp )
  outC = matrix( NA , nComp , nComp )
  out = matrix( NA , nComp , nComp )
  
  for( tel1 in 1:nComp )
  {
    for( tel2 in 1:nComp )
    {
      outA[ tel1 , tel2 ] = CosAngleVectors( A[,tel1] , A[,tel2] )
      outB[ tel1 , tel2 ] = CosAngleVectors( B[,tel1] , B[,tel2] )
      outC[ tel1 , tel2 ] = CosAngleVectors( C[,tel1] , C[,tel2] )
    }
  }
  tempM = outA * outB * outC
  tempEIG = eigen( tempM )
  
  out = list()
  out$matrix = tempM
  out$crit1 = sum( tempM <= -.85 ) # number of triple cosine products smaller than -.85 (is problematic)
  out$crit1_ok = sum( tempM <= -.85 ) == 0
  out$crit2 = min( tempEIG$values ) #smallest eigenvalue smaller than .50 (is problematic)
  out$crit2_ok = min( tempEIG$values ) >= .50
  out$crit3 = max( tempEIG$values ) / min( tempEIG$values ) #condition number > 5 (is problematic)
  out$crit3_ok = ( max( tempEIG$values ) / min( tempEIG$values ) ) <= 5
  out$tripcosprod_ok = out$crit1_ok & out$crit2_ok & out$crit3_ok 
  return(out)
}