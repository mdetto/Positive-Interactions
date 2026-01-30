##R Code for computing IAAR and ISAR excluding big trees

LargeTreeExclusion=function(main,tr){ 
  
  #parameters setting
  allsp = sort(unique(main$sp))   #------------------------------------------------------step 1
  n = length(allsp)      # number of species
  m = length(main$dist) # number of distanced analyzed
  dmax = 60         # maximum distance (m) for null model
  
 
  ## group individuals according to 20 dbh classes
  class=seq(0,log10(max(main$dbh)),length.out=21)
  group=numeric(length(main$gx))
  for (i in 1:20){
    group[log10(main$dbh)>=class[i] & log10(main$dbh)<class[i+1]]=i
  }
  
  ## pre-allocate outputs
  RNN = RNS = matrix(data=NA,nrow=n,ncol=m)
  z.N = z.S = matrix(data=NA,nrow=n,ncol=m)
  p.N = p.S = matrix(data=NA,nrow=n,ncol=m)
  abund = numeric(n)
  
  ## MAIN loop over species--------------------------------------------------------step 2
  for (i in 1:n){
    cat(n-i,"\r")
    
    use = main$sp==allsp[i] & main$dbh<tr
    n1 = sum(use)
    abund[i]=sum(main$sp==allsp[i])
    
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
      
      count = 0
      ## loop across individuals of species i
      for (j in 1:n1){
        
  		 # non-focal species in the same size class
  		 use = main$sp!=allsp[i] & group==g1[j]
  		 x0 = main$gx[use]
  		 y0 = main$gy[use]
  		 N0 = main$N[,use]
  		 S0 = main$S[,use]
	  
        # select individuals with similar DBH-------------------------------step 3.1 select individuals within similar DBH
        dx = x0 - x1[j]
        dy = y0 - y1[j]
        d2 = dx*dx + dy*dy                                         # compute distance from focal individual
        use1 = which(d2<dmax^2 & N0[1,]>=0)                        # select individuals within 60 m radius and in the same dbh class of focal individual
       #use1 = which(abs(ld1-ld0[j])<db & N1[1,]>=0)               # un-comment for excluding the 60 m radius restriction
        
        
        if (length(use1)>1){
          IAAR0 = IAAR0 + rowMeans(N0[,use1], na.rm = TRUE)  # IAAR0 include all heterospecifcs in the 60m area   # Eq. 3(a)
          ISAR0 = ISAR0 + rowMeans(S0[,use1], na.rm = TRUE)  # ISAR0 include all heterospecifcs in the 60m area   # Eq. 3(b)
          
          count = count+1
          
        }else if (length(use1)==1){ # if only one neighbor found, no random permutation
          IAAR0 = IAAR0 + N0[,use1]
          ISAR0 = ISAR0 + S0[,use1]
          
          count = count+1
        }
        
      }
      
      IAAR0=IAAR0/count
      ISAR0=ISAR0/count
      
      
      # compute statistics               ------------------------------------------step 4
      
      # mean relative IAAR and ISAR
      RNN[i,] = IAAR/IAAR0 # Eq. 2(a)
      RNS[i,] = ISAR/ISAR0 # Eq. 2(b)
      
      
    }
  }
  
  # return results              ------------------------------------------step 5
  
  RNN.pN = cbind(allsp, RNN, abund)
  RNS.pS = cbind(allsp, RNS, abund)
  
  
  return(list(RNN.pN=RNN.pN,RNS.pS=RNS.pS))
  
}





