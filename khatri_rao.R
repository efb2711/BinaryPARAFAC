khatri_rao <- function(A, B)
{
  R <- ncol(A)
  kr <- matrix(0, nrow = nrow(A) * nrow(B), ncol = R)
  for (r in 1:R)
  {
    kr[, r] <- as.vector(outer(A[, r], B[, r]))
  }
  return(kr)
}