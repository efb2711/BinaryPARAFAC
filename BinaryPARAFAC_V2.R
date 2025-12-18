BinaryPARAFAC_V2 <- function( X , I , J , K , R , tolerance = 1e-06 ,
                              lambda_grid = 0 , 
                              num_max_iters=100 , OptimMethod="CG" ,
                              StartSeed=123456 , nMultiStart=10 , 
                              include_offsets = c(TRUE, TRUE, TRUE) )
{
  #X is 2D matrix
  # Variables para multi-start
  if( length(lambda_grid) > 1 )
  {
    print( "ERROR: It is better to only specify a single lambda value" )
  }
  Info = list()
  Info$Sizes = c(I,J,K)
  Info$Rank = R
  Info$tolerance = tolerance
  Info$lambda_grid = lambda_grid
  Info$num_max_iters = num_max_iters
  Info$OptimMethod = OptimMethod
  Info$StartSeed = StartSeed
  Info$nMultiStart = nMultiStart
  Info$include_offsets = include_offsets
  
  LossMultiStarts <- c()
  TimeMultiStarts <- c()
  IterMultiStarts <- c()
  LossVectorMultistarts = list()
  
  best_loss <- Inf
  best_lambda <- NULL
  best_params <- NULL
  best_lossvector <- NULL
  seeds = GenerateSeed(nMultiStart, StartSeed)
  
  for (lambda in lambda_grid)
  {
#    cat("\nProbando lambda =", lambda, "\n")
    
    for (start in 1:nMultiStart)
    {
      set.seed(seeds[start])
###why using sd=.001????? it gives start values close to 0 for all component matrices
      A = matrix(rnorm(I*R, sd = 1), I, R)
      B = matrix(rnorm(J*R, sd = 1), J, R)
      C = matrix(rnorm(K*R, sd = 1), K, R)
      
      # Initialize the three offsets of a single index (mu_a, mu_b, mu_c)
      mu_a = matrix(0, I, 1) # Offset vector for mode I
      mu_b = matrix(0, J, 1) # Offset vector for mode J
      mu_c = matrix(0, K, 1) # Offset vector for mode K
      
      lossvector <- c()
      Lnew <- Loss(A, B, C, mu_a, mu_b, mu_c, X, lambda, include_offsets)
      
#      cat("Iteración inicial: Loss =", Lnew, " ||A|| =", norm(A,"F"), " ||B|| =", norm(B,"F"), " ||C|| =", norm(C,"F"), "\n")
#      cat( c("Multistart",start,"/",nMultiStart, "\n") )
            
      lossvector <- c(lossvector, Lnew)
      iter <- 0
      err <- 1
      start.time <- Sys.time()
      
      while (err > tolerance && iter < num_max_iters) 
      {
        iter <- iter + 1
        Lold <- Lnew
        
#        cat( paste0("Multistart ",start," / ",nMultiStart," -- ","iteration ",iter," / ",num_max_iters, "\n") )
        
        # Optimizar A
        resA <- optim(par = as.vector(A),
                      fn = function(par) Loss(matrix(par, I, R), B, C, mu_a, mu_b, mu_c, X, lambda, include_offsets),
                      gr = function(par) as.vector(grLogRegARec(matrix(par, I, R), B, C, mu_a, mu_b, mu_c, X, lambda, include_offsets)),
                      method = OptimMethod)
        A <- matrix(resA$par, I, R)
        lossvector <- c(lossvector, resA$value)
#        cat("Iter", iter, "- A: Loss =", resA$value, "Convergencia =", resA$convergence, "\n")
        
        
        # Optimizar B
        resB <- optim(par = as.vector(B),
                      fn = function(par) Loss(A, matrix(par, J, R), C, mu_a, mu_b, mu_c, X, lambda, include_offsets), # <--- A, [matrix(par, J, R)], C
                      gr = function(par) as.vector(grLogRegBRec(A, matrix(par, J, R), C, mu_a, mu_b, mu_c, X, lambda, include_offsets)), # <--- Función y posición corregidas
                      method = OptimMethod)
        B <- matrix(resB$par, J, R)
        lossvector <- c(lossvector, resB$value)
#        cat("Iter", iter, "- B: Loss =", resB$value, "Convergencia =", resB$convergence, "\n")
        
        # Optimizar C
        resC <- optim(par = as.vector(C),
                      fn = function(par) Loss(A, B, matrix(par, K, R), mu_a, mu_b, mu_c, X, lambda, include_offsets), # <--- A, B, [matrix(par, K, R)]
                      gr = function(par) as.vector(grLogRegCRec(A, B, matrix(par, K, R), mu_a, mu_b, mu_c, X, lambda, include_offsets)), # <--- Función y posición corregidas
                      method = OptimMethod)
        C <- matrix(resC$par, K, R)
        lossvector <- c(lossvector, resC$value,-9999)
#        cat("Iter", iter, "- C: Loss =", resC$value, "Convergencia =", resC$convergence, "\n")
        
        # Optimize mu_a
        if (include_offsets[1])
        {
          res_mua <- optim(par = as.vector(mu_a),
                           fn = function(par) Loss(A, B, C, matrix(par, I, 1), mu_b, mu_c, X, lambda, include_offsets),
                           gr = function(par) as.vector(grLogRegMua(A, B, C, matrix(par, I, 1), mu_b, mu_c, X, lambda, include_offsets)),
                           method = OptimMethod)
          mu_a <- matrix(res_mua$par, I, 1)
          lossvector <- c(lossvector, res_mua$value)
#          cat("Iter", iter, "- mu_a: Loss =", res_mua$value, "Convergencia =", res_mua$convergence, "\n")
        }
        
        # Optimize mu_b
        if (include_offsets[2])
        {
          res_mub <- optim(par = as.vector(mu_b),
                           fn = function(par) Loss(A, B, C, mu_a, matrix(par, J, 1), mu_c, X, lambda, include_offsets),
                           gr = function(par) as.vector(grLogRegMub(A, B, C, mu_a, matrix(par, J, 1), mu_c, X, lambda, include_offsets)),
                           method = OptimMethod)
          mu_b <- matrix(res_mub$par, J, 1)
          lossvector <- c(lossvector, res_mub$value)
#          cat("Iter", iter, "- mu_b: Loss =", res_mub$value, "Convergencia =", res_mub$convergence, "\n")
        }
        
        # Optimize mu_c
        if (include_offsets[3])
        {
          res_muc <- optim(par = as.vector(mu_c),
                           fn = function(par) Loss(A, B, C, mu_a, mu_b, matrix(par, K, 1), X, lambda, include_offsets),
                           gr = function(par) as.vector(grLogRegMuc(A, B, C, mu_a, mu_b, matrix(par, K, 1), X, lambda, include_offsets)),
                           method = OptimMethod)
          mu_c <- matrix(res_muc$par, K, 1)
          lossvector <- c(lossvector, res_muc$value)
#          cat("Iter", iter, "- mu_c: Loss =", res_muc$value, "Convergencia =", res_muc$convergence, "\n")
        }
        
        if( any(include_offsets) )
        {
          lossvector <- c(lossvector,-8888)
        }
        
        Lnew <- Loss(A, B, C, mu_a, mu_b, mu_c, X, lambda, include_offsets)
        err <- abs(Lold - Lnew) / abs(Lold)
#        cat("Iter", iter, "Loss =", Lnew, " Error Relativo =", err, "\n")
        
      }
      
      end.time <- Sys.time()
      TimeUsed <- round(end.time - start.time, 2)
      
      # Guardar resultados multi-start
      LossMultiStarts <- c(LossMultiStarts, Lnew)
      TimeMultiStarts <- c(TimeMultiStarts, TimeUsed)
      IterMultiStarts <- c(IterMultiStarts, iter)
      LossVectorMultistarts[[start]] = lossvector
      
      if (Lnew < best_loss)
      {
        best_loss <- Lnew
        best_lambda <- lambda
        best_loss_lambdazero <- Loss(A, B, C, mu_a, mu_b, mu_c, X, lambda=0, include_offsets)
        best_params <- list(A = A, B = B, C = C, mu_a = mu_a, mu_b = mu_b, mu_c = mu_c)
        best_lossvector <- lossvector
      }
    }
  }
  
  tripcosprod = ComputeTripleCosineProduct( best_params$A , best_params$B , 
                                            best_params$C )
  
  # Compute AIC y BIC
  num_offsets_included = sum( include_offsets * c(I,J,K) )
  nEffPar <- (I + J + K - 2) * R + num_offsets_included 
  AIC <- 2 * best_loss + 2 * nEffPar
  BIC <- 2 * best_loss + log(I*J*K) * nEffPar
  
#  cat("\nMejor lambda encontrado:", best_lambda, "\n")
#  cat("AIC:", AIC, " BIC:", BIC, "\n")
  
  return(list(A = best_params$A,
              B = best_params$B,
              C = best_params$C,
              mu_a = best_params$mu_a,
              mu_b = best_params$mu_b,
              mu_c = best_params$mu_c,
              tripcosprod = tripcosprod,
              best_lambda = best_lambda,
              best_loss = best_loss,
              best_lossvector = best_lossvector,
              best_loss_lambdazero = best_loss_lambdazero,
              LossMultiStarts = LossMultiStarts,
              TimeMultiStarts = TimeMultiStarts,
              LossVectorMultistarts = LossVectorMultistarts,
              AIC = AIC,
              BIC = BIC,
              nEffPar = nEffPar,
              TotalSize = I*J*K,
              Info = Info))
}