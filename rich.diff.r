#program 1
#calculate the local richness difference.
#different dbh classes are calcultated here: tr=1, tr=5, tr=10, tr=20, respectively

rich.diff=function(tr,S,dist){          
  
  plot(0,0,type="l",col = "red",xlim=c(2,16),ylim=c(-0.5,0.5))
  abline(h=0,cex=2)
  
  #dist=seq(2,16,2)
  dist1=dist2=dist3=dist4=dist5=dist6=dist7=dist8=numeric()
  dist1.max=dist1.min=numeric()
  
  R=30
  use2=which(main$dbh>tr)
  n=length(use2)
  
  for (i in 1:sp.no){            
    cat(sp.no-i,"\r")
    
    oth=which(main$sp!=allsp[i] & main$dbh>tr)
    no=length(oth)
    
    use=which(main$sp==allsp[i] & main$dbh>tr)
    nl=length(use)
    
    y1=numeric()
    y2 <- array(0, dim=c(8,30))
    dummy=(colMeans(S[use,])-colMeans(S[oth,]))/dist^2/pi
    for (j in 1:R){
      rp=sample(n)
      y1=colMeans(S[use2[rp[1:nl]],])/dist^2/pi-colMeans(S[use2[rp[(nl+1):(nl+no)]],])/dist^2/pi
      y2[,j]=y1
    }
    
    
    points(dist,dummy-rowSums(y2)/R,type="l",col=rainbow(sp.no)[i],cex=2)
    dist1[i]=mean(as.matrix(dummy-rowSums(y2)/R)[1])
    dist2[i]=mean(as.matrix(dummy-rowSums(y2)/R)[1:2])
    dist3[i]=mean(as.matrix(dummy-rowSums(y2)/R)[1:3])
    dist4[i]=mean(as.matrix(dummy-rowSums(y2)/R)[1:4])
    dist5[i]=mean(as.matrix(dummy-rowSums(y2)/R)[1:5])
    dist6[i]=mean(as.matrix(dummy-rowSums(y2)/R)[1:6])
    dist7[i]=mean(as.matrix(dummy-rowSums(y2)/R)[1:7])
    dist8[i]=mean(as.matrix(dummy-rowSums(y2)/R)[1:8])
    

    dist1.max[i]=max(as.matrix(dummy-rowSums(y2)/R)[1])
    dist1.min[i]=min(as.matrix(dummy-rowSums(y2)/R)[1])
  }

  return(list(dist1.max=dist1.max,dist1.min=dist1.min,dist1=dist1,dist2=dist2,dist3=dist3,dist4=dist4,dist5=dist5,dist6=dist6,dist7=dist7,dist8=dist8))
}


##program 2
# Randomization with 60 m from focal individuals, for all species
#library(emdbook)

rich.diff.60m=function(tr,S,dist,max.dist){   

plot(0,0,type="l",col = "red",xlim=c(2,16),ylim=c(0.9,1.1))
abline(h=1,cex=2)
  
# divide individuals in dbh classes
group=numeric(length(main$gx))
class=lseq(1,max(main$dbh)+1,21)

for (i in 1:20)
  group[main$dbh>=class[i] & main$dbh<class[i+1]]=i

# main loop  
Y=matrix(0,ncol=sp.no,nrow=length(dist))
p=matrix(0,ncol=sp.no,nrow=length(dist))
for (i in 1:sp.no){
  cat(sp.no-i,"\r")
  
  use=which(main$sp!=allsp[i]) #find all non-focal individuals
  x2=main$gx[use]
  y2=main$gy[use]
    
  use1=which(main$dbh>=tr & main$sp==allsp[i])
  if (length(use1)>1){  #if less than 2 stems, no need t.test
  x1=main$gx[use1]
  y1=main$gy[use1]
  H=matrix(0,nrow=length(use1),ncol=length(dist))
  for (j in 1:length(use1)){
    
    use2=which(group[use]==group[use1[j]])
    r2=(x2[use2]-x1[j])^2+(y2[use2]-y1[j])^2
    H[j,]=colMeans(S[use[use2[r2<max.dist^2]],])   
  }
  # for large dbh, sometime no "r2>max.dist^2", need to improve max.dist or exclude some individuals
  H=H[!rowSums(!is.finite(H)),] # H with NA values excluded, another case: note "m[!is.finite(m)] <- 0"
  
  Y[,i]=colMeans(S[use1,])/colMeans(as.data.frame(H))
  
  for (k in 1:length(dist))
    p[k,i]=t.test(S[use1,k],H[,k])$p.value   #even there is NA values, t.test no problem
  
  }
  
  points(dist,Y[,i],type="l",col=rainbow(sp.no)[i],cex=2)
  
 }
  return(list(Y=Y,p=p))
}



