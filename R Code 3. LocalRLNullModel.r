# R Code 4 | Script for local random labeling test, “LocalRLNullModel”.
# This null model takes into account the size structure of the community (by randomizing trees within the same size class), 
# and the habitat heterogeneity (by comparing individuals within a 60 m distance).

NullModel=function(main){ 
  
  #parameters setting
  allsp = sort(unique(main$sp))   #------------------------------------------------------step 1
  n = length(allsp)      # number of species
  m = length(main$dist) # number of distanced analyzed
  dmax = 60         # maximum distance (m) for null model
  
  ## group individuals according to 20 dbh classes
  class=seq(0,log10(max(main$dbh)),length.out=21)
  group=numeric(length(main$gx))
  for (i in 1:20)
    group[log10(main$dbh)>=class[i] & log10(main$dbh)<class[i+1]]=i
  
  group[log10(main$dbh)==log10(max(main$dbh))]=20
  
  ## pre-allocate outputs
  RNN = RNS = matrix(data=NA,nrow=n,ncol=m)
  z.N = z.S = matrix(data=NA,nrow=n,ncol=m)
  p.N = p.S = matrix(data=NA,nrow=n,ncol=m)
  abund = numeric(n)
  
  ## MAIN loop over species--------------------------------------------------------step 2
  for (i in 1:n){
    cat(n-i,"\r")
    
    use = main$sp==allsp[i]
    n1 = sum(use)
    abund[i]=n1
    
    if (n1>1){        # only consider species with more than 1 individual
      
      # focal species
      x1 = main$gx[use]
      y1 = main$gy[use]
      g1 = group[use]
      N1 = main$N[,use]
      S1 = main$S[,use]
      
      
      # compute ISAR---------------------------------------------------------------------step 3
      IAAR = rowMeans(N1, na.rm = TRUE)  # Eq. 1(a)
      ISAR = rowMeans(S1, na.rm = TRUE)  # Eq. 1(b)
      
      IAAR0 = ISAR0 = numeric(m)
      IAAR0.null = ISAR0.null = matrix(data=0,nrow=m,ncol=999)
      
      count = 0
      ## loop across individuals of species i
      for (j in 1:n1){
        
        # non-focal species within the same dbh class
        use = main$sp!=allsp[i] & group==g1[j]
        x0 = main$gx[use]
        y0 = main$gy[use]
        g0 = group[use]
        N0 = main$N[,use]
        S0 = main$S[,use]
        
        # select individuals with similar DBH-------------------------------step 3.1 select individuals within similar DBH
        dx = x0 - x1[j]
        dy = y0 - y1[j]
        d2 = dx*dx + dy*dy                                         # compute distance from focal individual
        use0 = which(d2<dmax^2 & g0==g1[j] & N0[1,]>=0)            # select individuals within 60 m radius and in the same dbh class of focal individual
        #use0 = which(g0==g1[j] & N0[1,]>=0)                        # un-comment for excluding the 60 m radius restriction
        
        
        if (length(use0)>1){
          IAAR0 = IAAR0 + rowMeans(N0[,use0], na.rm = TRUE)  # IAAR0 include all heterospecifcs in the 60m area   # Eq. 3(a)
          ISAR0 = ISAR0 + rowMeans(S0[,use0], na.rm = TRUE)  # ISAR0 include all heterospecifcs in the 60m area   # Eq. 3(b)
          
          k0 = sample(use0,999,replace = TRUE)   #------------------step 3.2 randomization (see methods: Test differences in neighborhood diversity by a random null model)
          IAAR0.null=IAAR0.null + N0[,k0]
          ISAR0.null=ISAR0.null + S0[,k0]
          count = count+1
          
        }else if (length(use0)==1){ # if only one neighbor found, no random permutation
          IAAR0 = IAAR0 + N0[,use0]
          ISAR0 = ISAR0 + S0[,use0]
          
          IAAR0.null=IAAR0.null+N0[,rep(use0,999)]
          ISAR0.null=ISAR0.null+S0[,rep(use0,999)]
          count = count+1
        }
        
      }
      
      IAAR0=IAAR0/count
      ISAR0=ISAR0/count
      
      IAAR0.null=IAAR0.null/count
      ISAR0.null=ISAR0.null/count
      
      # compute statistics               ------------------------------------------step 4
      
      # mean relative IAAR and ISAR
      RNN[i,] = IAAR/IAAR0 # Eq. 2(a)
      RNS[i,] = ISAR/ISAR0 # Eq. 2(b)
      
      # z-score
      
      z.N[i,] = (IAAR-rowMeans(IAAR0.null, na.rm = TRUE))/apply(IAAR0.null, 1, sd, na.rm = TRUE)
      z.S[i,] = (ISAR-rowMeans(ISAR0.null, na.rm = TRUE))/apply(ISAR0.null, 1, sd, na.rm = TRUE)
      
      
      # p-value
      for (k in 1:m) {
        # Calculate the two one-sided tails (including the observed x)
        # We add +1 to the numerator and denominator to account for the observed value
        p_low  <- (sum(IAAR0.null[k,] <= IAAR[k]) + 1) / 1000
        p_high <- (sum(IAAR0.null[k,] >= IAAR[k]) + 1) / 1000
        p.N[i,k] <- min(1.0, 2 * min(p_low, p_high))
        
        p_low  <- (sum(ISAR0.null[k,] <= ISAR[k]) + 1) / 1000
        p_high <- (sum(ISAR0.null[k,] >= ISAR[k]) + 1) / 1000
        p.S[i,k] <- min(1.0, 2 * min(p_low, p_high))
        
      }
    }
  }
  
  # return results              ------------------------------------------step 5
  
  RNN.pN = cbind(allsp, RNN, p.N, z.N, abund)
  RNS.pS = cbind(allsp, RNS, p.S, z.S, abund)
  
  
  return(list(RNN.pN=RNN.pN, RNS.pS=RNS.pS))
  
}





