#R Code 1 | Script for computing number of individuals (N) and number of species richness (S) around all stems, “ISAR”.
#compute number of individuals (N), number of species richness (S) of all individuals

# x,y:  stem coordinates
# Lx, Ly, size of the domain 
# Rmax is the max distance to focal stem that we want to calculate
# dr radial increment
# sp: species name
# dbh: diameter

ISAR <- function(x,y,sp,Lx,Ly,Rmax,dr){
  

  n<-length(x)
  R<-seq(dr,Rmax,dr)
  sp.num <- as.numeric(factor(sp))
  
  n.R = length(R)
  R2 = R*R
  Rmax2 = Rmax*Rmax
  
  # initialize
  S<-N<-matrix(0,n.R,n)
  
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
    
    
    dx = xi - xi[i]
    dy = yi - yi[i]
    r20 = dx*dx + dy*dy 
    
    use <- r20>0 & r20<=Rmax2
    z <- sum(use)
    if (z>0) {
      
      rj  <- r20[use]
      spj <- sp.num[use]
      
      S[n.R,i] <- length(unique(spj))
      N[n.R,i] <- z
      
    
    for (j in 1:n.R-1){
      
   
      use <- rj<=R2[j]
      z <- sum(use)
      
      if (z>0) {
        
        S[j,i] <- length(unique(spj[use]))
        N[j,i] <- z
        
      }
    }
    }
  }
  
 
  
  return(list(S = S,N = N, dist = R, n = no))
}