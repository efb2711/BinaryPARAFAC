is_binary_matrix <- function(matrix)
{
  return(all(matrix == 0 | matrix == 1))
}

#-------------------------------------------------------------------------------

sigmoide<-function(z)
{
  ( 1 / ( 1 + exp(-1*z) ) )
}

#-------------------------------------------------------------------------------

expand_offsets <- function(mu_a, mu_b, mu_c, I, J, K, include_offsets)
{
  contrib_a = 0
  if (include_offsets[1])
  {
    contrib_a <- mu_a %*% matrix(1, 1, J*K) # I x JK
  }
  
  contrib_b = 0
  if (include_offsets[2]) 
  {
    contrib_b_temp <- kronecker(matrix(1, K, 1), mu_b) # J*K x 1 
    contrib_b <- matrix(1, I, 1) %*% t(contrib_b_temp) # I x J*K
  }
  
  contrib_c = 0
  if (include_offsets[3])
  {
    contrib_c_temp <- kronecker(mu_c, matrix(1, J, 1)) # J*K x 1 
    contrib_c <- matrix(1, I, 1) %*% t(contrib_c_temp) # I x J*K
  }
  
  return(contrib_a + contrib_b + contrib_c)
}

#-------------------------------------------------------------------------------

Loss <- function(A, B, C, mu_a, mu_b, mu_c, Xa, lambda,include_offsets)
{
  I <- nrow(A)
  J <- nrow(B)
  K <- nrow(C)
  
  parafac_pred <- A %*% t(khatri_rao(C, B)) # I x JK
  mu_contrib <- expand_offsets(mu_a, mu_b, mu_c, I, J, K, include_offsets) 
  eta <- parafac_pred + mu_contrib
  Pia <- sigmoide(eta)
  
  eps <- 1e-12
  Pia <- pmin(pmax(Pia, eps), 1 - eps)
  L = sum(-1*Xa*log(Pia) - (1-Xa)*log(1-Pia))
  
  # Ridge penalty
  penalty_mua = ifelse(include_offsets[1], sum(mu_a^2), 0)
  penalty_mub = ifelse(include_offsets[2], sum(mu_b^2), 0)
  penalty_muc = ifelse(include_offsets[3], sum(mu_c^2), 0)
  penalty = lambda * (sum(A^2) + sum(B^2) + sum(C^2) + penalty_mua + penalty_mub + penalty_muc)
  
  return(L + penalty)
}

#-------------------------------------------------------------------------------

grLogRegARec <- function(A, B, C, mu_a, mu_b, mu_c, Xa, lambda, include_offsets)
{
  I <- nrow(A)
  J <- nrow(B)
  K <- nrow(C)
  
  parafac_pred <- (A %*% t(khatri_rao(C, B)))
  mu_contrib <- expand_offsets(mu_a, mu_b, mu_c, I, J, K,include_offsets)
  Pia <- sigmoide(parafac_pred + mu_contrib)
  
  gradA = (Pia - Xa) %*% khatri_rao(C, B) + 2 * lambda * A
  return(gradA)
}

#-------------------------------------------------------------------------------

grLogRegBRec <- function(A, B, C, mu_a, mu_b, mu_c, Xa, lambda, include_offsets)
{
  I <- nrow(A)
  J <- nrow(B)
  K <- nrow(C)  
  
  parafac_pred <- (A %*% t(khatri_rao(C, B)))
  mu_contrib <- expand_offsets(mu_a, mu_b, mu_c, I, J, K,include_offsets)
  Pia <- sigmoide(parafac_pred + mu_contrib)
  
  Pib <- permnew(Pia,I,J,K)
  Xb = permnew(Xa,I,J,K)
  
  gradB = (Pib - Xb) %*% khatri_rao(A, C) + 2 * lambda * B
  return(gradB)
} 

#-------------------------------------------------------------------------------

grLogRegCRec <- function(A, B, C, mu_a, mu_b, mu_c, Xa, lambda,include_offsets)
{
  I <- nrow(A)
  J <- nrow(B)
  K <- nrow(C)  
  
  parafac_pred <- (A %*% t(khatri_rao(C, B)))
  mu_contrib <- expand_offsets(mu_a, mu_b, mu_c, I, J, K,include_offsets)
  Pia <- sigmoide(parafac_pred + mu_contrib)
  
  Xb = permnew(Xa, I, J, K)
  Xc = permnew(Xb, J, K, I)
  Pib = permnew(Pia, I, J, K)
  Pic = permnew(Pib, J, K, I)
  
  gradC = (Pic - Xc) %*% khatri_rao(B, A) + 2 * lambda * C
  return(gradC)
}

#-------------------------------------------------------------------------------

grLogRegMua <- function(A, B, C, mu_a, mu_b, mu_c, Xa, lambda, include_offsets)
{
  I <- nrow(A)
  J <- nrow(B)
  K <- nrow(C)  
  
  parafac_pred <- (A %*% t(khatri_rao(C, B)))
  mu_contrib <- expand_offsets(mu_a, mu_b, mu_c, I, J, K,include_offsets)
  Pia <- sigmoide(parafac_pred + mu_contrib)
  
  grad = rowSums(Pia - Xa)
  
  penalty_term = 2 * lambda * mu_a
  
  return(matrix(grad, ncol = 1) + penalty_term)
}

#-------------------------------------------------------------------------------

grLogRegMub <- function(A, B, C, mu_a, mu_b, mu_c, Xa, lambda,include_offsets)
{
  I <- nrow(A)
  J <- nrow(B)
  K <- nrow(C)
  
  parafac_pred <- (A %*% t(khatri_rao(C, B)))
  mu_contrib <- expand_offsets(mu_a, mu_b, mu_c, I, J, K,include_offsets)
  Pia <- sigmoide(parafac_pred + mu_contrib)
  
  Xb = permnew(Xa, I, J, K) 
  Pib = permnew(Pia, I, J, K)
  
  grad = rowSums(Pib - Xb)
  
  penalty_term = 2 * lambda * mu_b
  
  return(matrix(grad, ncol = 1) + penalty_term)
}

#-------------------------------------------------------------------------------

grLogRegMuc <- function(A, B, C, mu_a, mu_b, mu_c, Xa, lambda,include_offsets)
{
  I <- nrow(A)
  J <- nrow(B)
  K <- nrow(C)
  
  parafac_pred <- (A %*% t(khatri_rao(C, B)))
  mu_contrib <- expand_offsets(mu_a, mu_b, mu_c, I, J, K,include_offsets)
  Pia <- sigmoide(parafac_pred + mu_contrib)
  
  Xb = permnew(Xa, I, J, K)
  Xc = permnew(Xb, J, K, I)
  Pib = permnew(Pia, I, J, K)
  Pic = permnew(Pib, J, K, I)
  
  grad = rowSums(Pic - Xc)
  
  penalty_term = 2 * lambda * mu_c
  
  return(matrix(grad, ncol = 1) + penalty_term)
}