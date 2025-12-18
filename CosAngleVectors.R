CosAngleVectors <- function( avec , bvec )
{
  # returns cosine of angle between two vectors avec and bvec
  #  (1 is perfect recovery, 0 no recovery, -1 is also possible)
  # avec (1 x nElem): vector
  # bvec (1 x nElem): vector
  div = sqrt( sum( avec^2 ) ) * sqrt( sum( bvec^2 ) )
  return( ( t(avec) %*% bvec ) / div )
}