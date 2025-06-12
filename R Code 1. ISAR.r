#R Code 1 | Script for computing number of individuals (N) and number of species richness (S) around all stems, “ISAR”.
#compute number of individuals (N), number of species richness (S), basal area (B) and shannon entropy (H) around all stems

# x,y:  stem coordinates
# Lx, Ly, size of the domain 
# Rmax is the max distance to focal stem that we want to calculate
# dr radial increment
# sp: species name
# dbh: diameter

ISAR <- function(x,y,sp,dbh,Lx,Ly,Rmax,dr){
  

  n<-length(x)
  R<-seq(dr,Rmax,dr)
  log.ba <- log(pi/4*dbh^2)
  mb = mean(log.ba)
  sp.num <- as.numeric(factor(sp))
  
  n.R = length(R)
  R2 = R*R
  Rmax2 = Rmax*Rmax
  
  # initialize
  S<-N<-B<-H<-matrix(0,n.R,n)
  
  # main loop across all individuals
  for (i in 1:n){
    cat(n-i,"\r")
    
  
    xi<-x
    yi<-y

    # periodic conditions
    if  (x[i]>Lx-Rmax){
      xi[xi<Rmax]<-xi[xi<Rmax]+Lx
    } else if (x[i]<Rmax){
      xi[xi>Lx-Rmax]<-xi[xi>Lx-Rmax]-Lx
    }

    if (y[i]<Rmax){
      yi[yi>Ly-Rmax]<-yi[yi>Ly-Rmax]-Ly
    } else if (y[i]>Ly-Rmax){
      yi[yi<Rmax]<-yi[yi<Rmax]+Ly
    }
    
    r20<-(xi-xi[i])^2+(yi-yi[i])^2
    
    use <- r20>0 & r20<=Rmax2
    z <- sum(use)
    if (z>0) {
      
      rj  <- r20[use]
      spj <- sp.num[use]
      log.baj <- log.ba[use]
      
      S[n.R,i] <- length(unique(spj))
      N[n.R,i] <- z
      B[n.R,i] <- sum(log.baj)/mb
      
      cc<-table(spj)/N[n.R,i]
      h<- -cc*log(cc)
      h[is.na(h)]<-0
      
      H[n.R,i] <- sum(h)   # Shannon entropy
      
    
    for (j in 1:n.R-1){
      
   
      use <- rj<=R2[j]
      z <- sum(use)
      
      if (z>0) {
        
        S[j,i] <- length(unique(spj[use]))
        N[j,i] <- z
        B[j,i] <- sum(log.baj[use])/mb
        
        cc<-table(spj[use])/N[j,i]
        h<- -cc*log(cc)
        h[is.na(h)]<-0
        
        H[j,i] <- sum(h)   # Shannon entropy
        
      }
    }
    }
  }
  
  
  species <- unique(sp.num)
  sp.no   <- length(species)
  
  SIR.S <- SIR.N <- SIR.B <- SIR.H <- matrix(0,length(R),sp.no)

  no <- numeric(sp.no)
  for (i in 1:sp.no) {
    
    use <- sp.num==species[i]
    no[i] <- sum(use)
    if (no[i] >1){
      SIR.S[,i] <- rowMeans(S[,use])
      SIR.B[,i] <- rowMeans(B[,use])
      SIR.N[,i] <- rowMeans(N[,use])
      SIR.H[,i] <- rowMeans(H[,use])
    } else if (no[i]==1){
      SIR.S[,i] <- S[,use]
      SIR.B[,i] <- B[,use]
      SIR.N[,i] <- N[,use]
      SIR.H[,i] <- H[,use]
    }
  }
  
  
  return(list(S = S,N = N,B = B,H = H,dist = R, n = no,
              SIR.S = SIR.S,
              SIR.N = SIR.N,
              SIR.B = SIR.B,
              SIR.H = SIR.H))
}