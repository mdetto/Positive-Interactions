#R Code 3 | Script for analytical tests based on two null models, Random Labeling (RL) model and Complete Spatial Randomness (CSR) model, “AnalyticalNullModel”.

NullModel <- function(n, E0.S, E0.N, E0.B, cvb, type, ...) {
  # Description
  # NullModel(res, cvb, 'CSR', lambda)
  #
  #   res:   output from ISAR.r
  #   cvb: log(basal area) coefficient of variation (std / avg) across all individuals
  #   lambda: intensity of Poisson point process (N/area * pi * radius^2 )
  #
  # NullModel(n, cvb, 'shuffle', Pk, k)
  #
  #   n:   vector of species abundance
  #   cvb: log(basal area) coefficient of variation (std / avg) across all individuals
  #   Pk:  probability of a neighborhood with k individuals (vector)
  #   k:   vector of k
  
  
  p <- n / sum(n)
  S <- length(n)
  
  if (type == "CSR") {
    # Complete Spatial Randomness
    lambda <- list(...)[[1]]

    C <- 0
    for (i in 1:S) {
      for (j in 1:S) {
        if (j != i) {
          C <- C + 1 + exp(-lambda * (p[i] + p[j])) - exp(-lambda * p[i]) - exp(-lambda * p[j])
        }
      }
    }
    
    E.S <- S - sum(exp(-lambda * p))
    V.S <- C + E.S * (1 - E.S)
    
    E.B <- lambda
    V.B <- lambda * (cvb^2 + 1)
    
    E.N <- lambda
    V.N <- lambda
    
  } else if (type == "shuffle") {
    # Species Shuffle
    Pk <- list(...)[[1]]
    k <- list(...)[[2]]
    K <- length(k)
    Ek <- numeric(K)
    Ck <- numeric(K)
    
    for (i in 1:K) {
      kk <- k[i]
      Ek[i] <- sum(1 - (1 - p)^kk)
      
      C <- 0
      for (m in 1:S) {
        for (l in 1:S) {
          if (l != m) {
            C <- C + 1 - (1 - p[m])^kk - (1 - p[l])^kk + (1 - p[m] - p[l])^kk
          }
        }
      }
      
      Ck[i] <- C
    }
    
    E.S <- sum(Ek * Pk)
    V.S <- sum(Ck * Pk) + E.S * (1 - E.S)
    
    mk <- sum(Pk * k)
    vk <- sum(Pk * k^2) - mk^2
    E.B <- mk
    V.B <- mk * cvb^2 + vk
    
    E.N <- mk
    V.N <- vk
  }
  
  # Confidence intervals for species richness
  alfa <- 0.95
  ci.S <- matrix(0, nrow = S, ncol = 2)
  ci.S[, 1] <- E.S - sqrt(V.S / n) * sqrt(2) * qnorm(alfa)
  ci.S[, 2] <- E.S + sqrt(V.S / n) * sqrt(2) * qnorm(alfa)
  
  # Confidence intervals for species abundnace
  ci.N <- matrix(0, nrow = S, ncol = 2)
  ci.N[, 1] <- E.N - sqrt(V.N / n) * sqrt(2) * qnorm(alfa)
  ci.N[, 2] <- E.N + sqrt(V.N / n) * sqrt(2) * qnorm(alfa)
  
  # Confidence intervals for species abundnace (basal area)
  ci.B <- matrix(0, nrow = S, ncol = 2)
  ci.B[, 1] <- E.B - sqrt(V.B / n) * sqrt(2) * qnorm(alfa)
  ci.B[, 2] <- E.B + sqrt(V.B / n) * sqrt(2) * qnorm(alfa)
  
  # two-tail test
  z.score = (E0.S-E.S)/sqrt(V.S/n)
  p_value.S <- 2 * pnorm(abs(z.score), lower.tail = FALSE)
  
  z.score = (E0.N-E.N)/sqrt(V.N/n)
  p_value.N <- 2 * pnorm(abs(z.score), lower.tail = FALSE)
  
  z.score = (E0.B-E.B)/sqrt(V.B/n)
  p_value.B <- 2 * pnorm(abs(z.score), lower.tail = FALSE)
  
  
  return(list(ci.S = ci.S, ci.N = ci.N, ci.B = ci.B, 
              p_value.S = p_value.S, p_value.N = p_value.N, p_value.B = p_value.B,
              E.S = E.S, E.N = E.N,  E.B = E.B, 
              V.S = V.S, V.N = V.N, V.B = V.B))
}
