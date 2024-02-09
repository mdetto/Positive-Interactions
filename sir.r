#compute number of individuals (N), number of species richness (S), basal area (B) and shannon equivalent richness (H)around all stems

# x,y  point positions
# dr, width of the anulii
# Lx, Ly, size of the domain 
# Rmax is the max distance to focal stem that we want to calculate

sir=function(x,y,sp,dbh,dr,Lx,Ly,Rmax){
  
  R=seq(0,Rmax,dr)
  n=length(x)

# % all individuals
  S=N=B=H=matrix(0,length(R)-1,n)

  for (i in 1:n){
	  cat(n-i,"\r")
    
	  for (j in 1:(length(R)-1)){
	
	# focal individual
	xi=x
        yi=y
	  # periodic conditions
	  if  (x[i]>Lx-R[j+1]){
		xi[xi<R[j+1]]=xi[xi<R[j+1]]+Lx
	  }

	  if (y[i]>Ly-R[j+1]){
	  yi[yi<R[j+1]]=yi[yi<R[j+1]]+Ly
	  }

	  if (x[i]<R[j+1]){
	  xi[xi>Lx-R[j+1]]=xi[xi>Lx-R[j+1]]-Lx
	  }

	  if (y[i]<R[j+1]){
	  yi[yi>Ly-R[j+1]]=yi[yi>Ly-R[j+1]]-Ly
   	}


	  r=sqrt((xi-xi[i])^2+(yi-yi[i])^2)
	  use=which(r>0 & r<=R[j+1])
	  if (length(use)>0) {
      
	  S[j,i]=length(unique(sp[use]))
          N[j,i]=length(use)
	  B[j,i]=sum(dbh[use]^2)
	  
          cc=table(sp[use])/N[j,i]
	  h= -cc*log(cc)
	  h[is.na(h)]<-0
	  H[j,i] = exp(sum(h))   # shannon equivalent richness
	  
	  }
	}
}

 B=pi/4*B
 dist=seq(2,Rmax,dr)
 return(list(S=S,N=N,B=B,H=H,dist=dist))
 }
