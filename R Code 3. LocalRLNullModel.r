##R Code 4 | Script for local random labelling test, “LocalRLNullModel”.
##This null model takes into account the size structure of the community (by randomizing trees within the same size class), and the habitat heterogeneity (by comparing individuals within a 60 m distance).

rich.diff.ran.sp1=function(mainr){ 

#parameters setting
allsp = sort(unique(mainr$sp))   #------------------------------------------------------step 1
n = length(allsp) # number of species
m = length(res$dist) # number of distanced analyzed

## parameters setting
dmax = 60        # maximum distance (m) for null model
db   = 0.1       # log10 DBH interval for size classes

################# randomization method 2 #########################
## MAIN loop over species--------------------------------------------------------step 2
RNN = RNS = RNB = RNH = matrix(data=NA,nrow=n,ncol=m)
RNN.dn = RNS.dn = RNB.dn = RNH.dn = matrix(data=NA,nrow=n,ncol=m)
RNN.up = RNS.up = RNB.up = RNH.up = matrix(data=NA,nrow=n,ncol=m)
p.N=p.S=p.B=p.H=matrix(0,ncol=n,nrow=m)

for (i in 1:n){
  cat(n-i,"\r")
 
  use = mainr$sp==allsp[i]
  n1 = sum(use)
 
  if (n1>1){        # only consider species with more than 1 individual

    # focal species
    x0 = mainr$gx[use]
    y0 = mainr$gy[use]
    ld0 = log10(mainr$dbh[use])
    S.N = res$N[,use]
    S.S = res$S[,use]
    S.B = res$B[,use]
    S.H = res$H[,use]

    # non-focal species
    x1 = mainr$gx[!use]
    y1 = mainr$gy[!use]
    ld1 = log10(mainr$dbh[!use])
    N1 = res$N[,!use]
    S1 = res$S[,!use]
    B1 = res$B[,!use]
    H1 = res$H[,!use]
   
    # compute ISAR---------------------------------------------------------------------step 2.1
    ISAR.N = rowMeans(S.N, na.rm = TRUE)
    ISAR.S = rowMeans(S.S, na.rm = TRUE)
    ISAR.B = rowMeans(S.B, na.rm = TRUE)
    ISAR.H = rowMeans(S.H, na.rm = TRUE)
    
    # compute ISAR.null----------------------------------------------------------------step 2.2   
    ISAR.null.N = ISAR.null.S = ISAR.null.B = ISAR.null.H = matrix(data=NA,nrow=m,ncol=999)
    H.N = H.S = H.B = H.H = array(data=NA, c(n1,m,999))

    for (j in 1:n1){
     
      # select individuals with similar DBH-------------------------------step 2.2 select individuals within similar DBH
      d = sqrt((x1-x0[j])^2+(y1-y0[j])^2)                      #compute distance from focal individual
      use1 = which(d<dmax & abs(ld1-ld0[j])<db)                #select individuals within 60 m radius and in the same dbh class of focal individual
      #use1 = which(abs(ld1-ld0[j])<db)                        #without considering 60 m radius restriction

      if (length(use1)>0) {     #-----------------------------------------this "0" could be changed to 1,2, or larger values if problems occur
       
        k0 = sample(use1,999,replace = TRUE)   #-------------------------------------step 2.2 species randomization
        H.N[j,,] = N1[,k0]
        H.S[j,,] = S1[,k0]
        H.B[j,,] = B1[,k0]
        H.H[j,,] = H1[,k0]
      }
    }
      
    ISAR.null.N = apply(H.N, 3, colMeans, na.rm = TRUE);
    ISAR.null.S = apply(H.S, 3, colMeans, na.rm = TRUE);
    ISAR.null.B = apply(H.B, 3, colMeans, na.rm = TRUE);
    ISAR.null.H = apply(H.H, 3, colMeans, na.rm = TRUE);

    # compute mean relative ISAR-----------------------------------------------------step 2.3
    RNN[i,] = ISAR.N/rowMeans(ISAR.null.N, na.rm = TRUE)
    RNS[i,] = ISAR.S/rowMeans(ISAR.null.S, na.rm = TRUE)
    RNB[i,] = ISAR.B/rowMeans(ISAR.null.B, na.rm = TRUE)
    RNH[i,] = ISAR.H/rowMeans(ISAR.null.H, na.rm = TRUE)   
      
    # take 0.025 and 0.975 quantiles and std-------------------------------------------step 2.4 
    RNN.dn[i,] = ISAR.N/apply(ISAR.null.N, 1, quantile, probs = .025, na.rm = TRUE )
    RNS.dn[i,] = ISAR.S/apply(ISAR.null.S, 1, quantile, probs = .025, na.rm = TRUE )
    RNB.dn[i,] = ISAR.B/apply(ISAR.null.B, 1, quantile, probs = .025, na.rm = TRUE )
    RNH.dn[i,] = ISAR.H/apply(ISAR.null.H, 1, quantile, probs = .025, na.rm = TRUE )
   
    RNN.up[i,] = ISAR.N/apply(ISAR.null.N, 1, quantile, probs = .975, na.rm = TRUE )
    RNS.up[i,] = ISAR.S/apply(ISAR.null.S, 1, quantile, probs = .975, na.rm = TRUE )
    RNB.up[i,] = ISAR.B/apply(ISAR.null.B, 1, quantile, probs = .975, na.rm = TRUE )
    RNH.up[i,] = ISAR.H/apply(ISAR.null.H, 1, quantile, probs = .975, na.rm = TRUE )

    RNN.sd[i,] = ISAR.N/apply(ISAR.null.N, 1, sd, probs = .975, na.rm = TRUE )
    RNS.sd[i,] = ISAR.S/apply(ISAR.null.S, 1, sd, probs = .975, na.rm = TRUE )
    RNB.sd[i,] = ISAR.B/apply(ISAR.null.B, 1, sd, probs = .975, na.rm = TRUE )
    RNH.sd[i,] = ISAR.H/apply(ISAR.null.H, 1, sd, probs = .975, na.rm = TRUE )

  }
    #add t.test----------------------------------------------------------------------step 2.5
    for (k in 1:m) {
    p.N[k,i]=t.test(S.N[k,],ISAR.null.N[k,])$p.value
    p.S[k,i]=t.test(S.S[k,],ISAR.null.S[k,])$p.value
    p.B[k,i]=t.test(S.B[k,],ISAR.null.B[k,])$p.value
    p.H[k,i]=t.test(S.H[k,],ISAR.null.H[k,])$p.value}
}

  RNN.pN=cbind(allsp, RNN,t(p.N),RNN.dn,RNN.up,RNN.sd)
  RNS.pS=cbind(allsp, RNS,t(p.S),RNS.dn,RNS.up,RNS.sd)
  RNB.pB=cbind(allsp, RNB,t(p.B),RNB.dn,RNB.up,RNB.sd)
  RNH.pH=cbind(allsp, RNB,t(p.H),RNB.dn,RNB.up,RNH.sd)
 
  return(list(RNN.pN=RNN.pN,RNS.pS=RNS.pS))

}





