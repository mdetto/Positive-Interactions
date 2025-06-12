R Code 2 | Script for calculating the relative neighborhood abundance and richness for trees within 60 m circular neighborhood of focal species and falling within the same diameter at breast height (DBH) class as the focal tree, “RNN.RNS”.
#library(emdbook) needed

rich.diff.60m=function(tr,S,dist){   

plot(0,0,type="l",col = "red",xlim=c(2,16),ylim=c(0.9,1.1))
abline(h=1,cex=2)
  
# divide individuals in dbh classes
group=numeric(length(main$gx))
class=lseq(1,max(main$dbh),21)

for (i in 1:20)
  group[main$dbh>=class[i] & main$dbh<class[i+1]]=i

# main loop  
Y=matrix(0,ncol=sp.no,nrow=length(dist))
p=matrix(0,ncol=sp.no,nrow=length(dist))
for (i in 1:sp.no){
  cat(sp.no-i,"\r")
  
  use=which(main$sp!=allsp[i])      #find all non-focal individuals
  x2=main$gx[use]
  y2=main$gy[use]
    
  use1=which(main$dbh>=tr & main$sp==allsp[i])

  if (length(use1)>start.n){   
  x1=main$gx[use1]
  y1=main$gy[use1]
  H=matrix(0,nrow=length(use1),ncol=length(dist))
  for (j in 1:length(use1)){
    
    use2=which(group[use]==group[use1[j]])
    r2=(x2[use2]-x1[j])^2+(y2[use2]-y1[j])^2
    H[j,]=colMeans(S[use[use2[r2<max.dist^2]],]) 

  }

  H=H[!rowSums(!is.finite(H)),] # H with NA values excluded, another case: note "m[!is.finite(m)] <- 0"
  
  Y[,i]=colMeans(S[use1,])/colMeans(H)
  
  for (k in 1:length(dist))
    p[k,i]=t.test(S[use1,k],H[,k])$p.value
  
  }
  
  points(dist,Y[,i],type="l",col=rainbow(sp.no)[i],cex=2)
  
 }
  return(list(Y=Y,p=p))
}


