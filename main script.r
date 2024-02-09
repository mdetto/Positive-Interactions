#1
#read data file and set the default folder
#assume that the big plot data file "data" contains four columns: gx, gy, sp and dbh.
#set the names of data,"data" indicates the plot data file.

setwd('d:/rich')  #where all result excel files exported and stored.

data=read.csv("d:/data.csv",header=T)   #need the original data file for the following neighbor richness calculation and test

names(data)[1]="sp"          # "1" can be changed by your column number of "sp"  in the file
names(data)[2]="gx"          # "4" can be changed by your column number of coordiante "gx"  in the file
names(data)[3]="gy"          # "5" can be changed by your column number of coordiante "gy"  in the file
names(data)[4]="dbh"         # "3" can be changed by your column number of "dbh"  in the file

main=data


#2
#define the sp names of all species

allsp=sort(unique(main$sp))
sp.no=length(allsp)

write.csv(allsp, "d:/rich/allsp.csv")   #species name list


#3 calculate the number of individuals and basal area for all species and legumes

No.stems1=tapply(main$dbh,main$sp,length) # 1cm
No.ba1=tapply(pi*main$dbh^2/4,main$sp,sum) #1cm
No.stems5=tapply(main$dbh[main$dbh>=50],main$sp[main$dbh>=50],length) # 5cm
No.ba5=tapply(pi*main$dbh[main$dbh>=50]^2/4,main$sp[main$dbh>=50],sum) #5cm
No.stems10=tapply(main$dbh[main$dbh>=100],main$sp[main$dbh>=100],length) # 10cm
No.ba10=tapply(pi*main$dbh[main$dbh>=100]^2/4,main$sp[main$dbh>=100],sum) #10cm
No.stems20=tapply(main$dbh[main$dbh>=200],main$sp[main$dbh>=200],length) # 20cm
No.ba20=tapply(pi*main$dbh[main$dbh>=200]^2/4,main$sp[main$dbh>=200],sum) #20cm

No.stems=cbind(No.stems1,No.stems5,No.stems10,No.stems20)
No.ba=cbind(No.ba1,No.ba5,No.ba10,No.ba20)

write.csv(No.stems, "d:/rich/No.stem.csv")  #no of stems for all species with four start DBH size
write.csv(No.ba, "d:/rich/No.ba.csv")       #total basal area for all species with four start DBH size


#4
#read the "sir" file to compute number of species richness (S), number of individuals (N), basal area (B) and shannon equivalent richness (H) around all focal stems

source('d:/rich/sir.r')

# res=sir(main$gx,main$gy,main$sp,main$dbh,2,20,20,16)  #by this, we get S, N, B, H
# S=t(res$S)
# N=t(res$N)
# B=t(res$B)
# H=t(res$H)
# write.csv(S,"d:/Plotname.S.csv")
# write.csv(N,"d:/Plotname.N.csv")
# write.csv(B,"d:/Plotname.B.csv")
# write.csv(H,"d:/Plotname.H.csv")
# S=as.data.frame(S)
# N=as.data.frame(N)
# B=as.data.frame(B)
# H=as.data.frame(H)

#if you have calculated the S, N, B, H last time, use the following codes to read S, N, B, H instead using the codes above in step 4.
S=read.csv("d:/rich/Plotname.S.csv",header=T)
S=as.data.frame(S)
S=S[,2:9]
names(S)

N=read.csv("d:/rich/Plotname.N.csv",header=T)
N=as.data.frame(N)
N=N[,2:9]
names(N)

B=read.csv("d:/rich/Plotname.B.csv",header=T)
B=as.data.frame(B)
B=B[,2:9]
names(B)

H=read.csv("d:/rich/Plotname.H.csv",header=T)
H=as.data.frame(H)
H=H[,2:9]
names(H)


#5
#including three parts: (1)all species by permutation, (2)all species with randomization
#calculate the local richness (S, N, B, H) difference.
#different dbh classes are calcultated here: tr=1 cm, tr=5 cm, tr=10 cm, tr=20 cm, respectively
#in case some species do not have large DBH sizes, so all results are saved seperately.
#those species without enough DBH sizes will not be calculated.

source('D:/rich/rich.diff.r')

install.packages("emdbook")    #need to install this package
library("emdbook")

max.dist=60   #max distance for randomization test, Radius around focal species
dist=seq(2,16,2)   #distance around the focal species from 2 to 16 m, with 2 m interval

##PART 1, rich.diff, permutation
##for S
rich.diff1=rich.diff(10,S,dist)  
dist.all1=cbind(as.data.frame(allsp),rich.diff1$dist1, rich.diff1$dist2, rich.diff1$dist3, rich.diff1$dist4, rich.diff1$dist5, rich.diff1$dist6, rich.diff1$dist7, rich.diff1$dist8)
write.csv(dist.all1,"d:/rich/dist.all1.S.csv")

rich.diff5=rich.diff(50,S,dist)  
dist.all5=cbind(as.data.frame(allsp),rich.diff5$dist1, rich.diff5$dist2, rich.diff5$dist3, rich.diff5$dist4, rich.diff5$dist5, rich.diff5$dist6, rich.diff5$dist7, rich.diff5$dist8)
write.csv(dist.all5,"d:/rich/dist.all5.S.csv")

rich.diff10=rich.diff(100,S,dist) 
dist.all10=cbind(as.data.frame(allsp),rich.diff10$dist1, rich.diff10$dist2, rich.diff10$dist3, rich.diff10$dist4, rich.diff10$dist5, rich.diff10$dist6, rich.diff10$dist7, rich.diff10$dist8)
write.csv(dist.all10,"d:/rich/dist.all10.S.csv")

rich.diff20=rich.diff(200,S,dist) 
dist.all20=cbind(as.data.frame(allsp),rich.diff20$dist1, rich.diff20$dist2, rich.diff20$dist3, rich.diff20$dist4, rich.diff20$dist5, rich.diff20$dist6, rich.diff20$dist7, rich.diff20$dist8)
write.csv(dist.all20,"d:/rich/dist.all20.S.csv")

##for N
rich.diff1=rich.diff(10,N,dist)  
dist.all1=cbind(as.data.frame(allsp),rich.diff1$dist1, rich.diff1$dist2, rich.diff1$dist3, rich.diff1$dist4, rich.diff1$dist5, rich.diff1$dist6, rich.diff1$dist7, rich.diff1$dist8)
write.csv(dist.all1,"d:/rich/dist.all1.N.csv")

rich.diff5=rich.diff(50,N,dist)  
dist.all5=cbind(as.data.frame(allsp),rich.diff5$dist1, rich.diff5$dist2, rich.diff5$dist3, rich.diff5$dist4, rich.diff5$dist5, rich.diff5$dist6, rich.diff5$dist7, rich.diff5$dist8)
write.csv(dist.all5,"d:/rich/dist.all5.N.csv")

rich.diff10=rich.diff(100,N,dist) 
dist.all10=cbind(as.data.frame(allsp),rich.diff10$dist1, rich.diff10$dist2, rich.diff10$dist3, rich.diff10$dist4, rich.diff10$dist5, rich.diff10$dist6, rich.diff10$dist7, rich.diff10$dist8)
write.csv(dist.all10,"d:/rich/dist.all10.N.csv")

rich.diff20=rich.diff(200,N,dist) 
dist.all20=cbind(as.data.frame(allsp),rich.diff20$dist1, rich.diff20$dist2, rich.diff20$dist3, rich.diff20$dist4, rich.diff20$dist5, rich.diff20$dist6, rich.diff20$dist7, rich.diff20$dist8)
write.csv(dist.all20,"d:/rich/dist.all20.N.csv")

##for B
rich.diff1=rich.diff(10,B,dist)  
dist.all1=cbind(as.data.frame(allsp),rich.diff1$dist1, rich.diff1$dist2, rich.diff1$dist3, rich.diff1$dist4, rich.diff1$dist5, rich.diff1$dist6, rich.diff1$dist7, rich.diff1$dist8)
write.csv(dist.all1,"d:/rich/dist.all1.B.csv")

rich.diff5=rich.diff(50,B,dist)  
dist.all5=cbind(as.data.frame(allsp),rich.diff5$dist1, rich.diff5$dist2, rich.diff5$dist3, rich.diff5$dist4, rich.diff5$dist5, rich.diff5$dist6, rich.diff5$dist7, rich.diff5$dist8)
write.csv(dist.all5,"d:/rich/dist.all5.B.csv")

rich.diff10=rich.diff(100,B,dist) 
dist.all10=cbind(as.data.frame(allsp),rich.diff10$dist1, rich.diff10$dist2, rich.diff10$dist3, rich.diff10$dist4, rich.diff10$dist5, rich.diff10$dist6, rich.diff10$dist7, rich.diff10$dist8)
write.csv(dist.all10,"d:/rich/dist.all10.B.csv")

rich.diff20=rich.diff(200,B,dist) 
dist.all20=cbind(as.data.frame(allsp),rich.diff20$dist1, rich.diff20$dist2, rich.diff20$dist3, rich.diff20$dist4, rich.diff20$dist5, rich.diff20$dist6, rich.diff20$dist7, rich.diff20$dist8)
write.csv(dist.all20,"d:/rich/dist.all20.B.csv")


##for H
rich.diff1=rich.diff(10,H,dist)  
dist.all1=cbind(as.data.frame(allsp),rich.diff1$dist1, rich.diff1$dist2, rich.diff1$dist3, rich.diff1$dist4, rich.diff1$dist5, rich.diff1$dist6, rich.diff1$dist7, rich.diff1$dist8)
write.csv(dist.all1,"d:/rich/dist.all1.H.csv")

rich.diff5=rich.diff(50,H,dist)  
dist.all5=cbind(as.data.frame(allsp),rich.diff5$dist1, rich.diff5$dist2, rich.diff5$dist3, rich.diff5$dist4, rich.diff5$dist5, rich.diff5$dist6, rich.diff5$dist7, rich.diff5$dist8)
write.csv(dist.all5,"d:/rich/dist.all5.H.csv")

rich.diff10=rich.diff(100,H,dist) 
dist.all10=cbind(as.data.frame(allsp),rich.diff10$dist1, rich.diff10$dist2, rich.diff10$dist3, rich.diff10$dist4, rich.diff10$dist5, rich.diff10$dist6, rich.diff10$dist7, rich.diff10$dist8)
write.csv(dist.all10,"d:/rich/dist.all10.H.csv")

rich.diff20=rich.diff(200,H,dist) 
dist.all20=cbind(as.data.frame(allsp),rich.diff20$dist1, rich.diff20$dist2, rich.diff20$dist3, rich.diff20$dist4, rich.diff20$dist5, rich.diff20$dist6, rich.diff20$dist7, rich.diff20$dist8)
write.csv(dist.all20,"d:/rich/dist.all20.H.csv")



##PART 2, rich.diff.60m, Y and p(t.test), all species
##for S
rich.diff1=rich.diff.60m(10,S,dist,max.dist)  
dist.all1.Yp=cbind(as.data.frame(allsp),t(rich.diff1$Y),t(rich.diff1$p))
write.csv(dist.all1.Yp,"d:/rich/dist.all1.S.Yp.csv")

rich.diff5=rich.diff.60m(50,S,dist,max.dist)  
dist.all5.Yp=cbind(as.data.frame(allsp),t(rich.diff5$Y),t(rich.diff5$p))
write.csv(dist.all5.Yp,"d:/rich/dist.all5.S.Yp.csv")

rich.diff10=rich.diff.60m(100,S,dist,max.dist) 
dist.all10.Yp=cbind(as.data.frame(allsp),t(rich.diff10$Y),t(rich.diff10$p))
write.csv(dist.all10.Yp,"d:/rich/dist.all10.S.Yp.csv")

rich.diff20=rich.diff.60m(200,S,dist,max.dist) 
dist.all20.Yp=cbind(as.data.frame(allsp),t(rich.diff20$Y),t(rich.diff20$p))
write.csv(dist.all20.Yp,"d:/rich/dist.all20.S.Yp.csv")

##for N
rich.diff1=rich.diff.60m(10,N,dist,max.dist)  
dist.all1.Yp=cbind(as.data.frame(allsp),t(rich.diff1$Y),t(rich.diff1$p))
write.csv(dist.all1.Yp,"d:/rich/dist.all1.N.Yp.csv")

rich.diff5=rich.diff.60m(50,N,dist,max.dist)  
dist.all5.Yp=cbind(as.data.frame(allsp),t(rich.diff5$Y),t(rich.diff5$p))
write.csv(dist.all5.Yp,"d:/rich/dist.all5.N.Yp.csv")

rich.diff10=rich.diff.60m(100,N,dist,max.dist) 
dist.all10.Yp=cbind(as.data.frame(allsp),t(rich.diff10$Y),t(rich.diff10$p))
write.csv(dist.all10.Yp,"d:/rich/dist.all10.N.Yp.csv")

rich.diff20=rich.diff.60m(200,N,dist,max.dist) 
dist.all20.Y=cbind(as.data.frame(allsp),t(rich.diff20$Y),t(rich.diff20$p))
write.csv(dist.all20.Yp,"d:/rich/dist.all20.N.Yp.csv")

##for B
rich.diff1=rich.diff.60m(10,B,dist,max.dist)  
dist.all1.Yp=cbind(as.data.frame(allsp),t(rich.diff1$Y),t(rich.diff1$p))
write.csv(dist.all1.Yp,"d:/rich/dist.all1.B.Yp.csv")

rich.diff5=rich.diff.60m(50,B,dist,max.dist)  
dist.all5.Yp=cbind(as.data.frame(allsp),t(rich.diff5$Y),t(rich.diff5$p))
write.csv(dist.all5.Yp,"d:/rich/dist.all5.B.Yp.csv")

rich.diff10=rich.diff.60m(100,B,dist,max.dist) 
dist.all10.Yp=cbind(as.data.frame(allsp),t(rich.diff10$Y),t(rich.diff10$p))
write.csv(dist.all10.Yp,"d:/rich/dist.all10.B.Yp.csv")

rich.diff20=rich.diff.60m(200,B,dist,max.dist) 
dist.all20.Yp=cbind(as.data.frame(allsp),t(rich.diff20$Y),t(rich.diff20$p))
write.csv(dist.all20.Yp,"d:/rich/dist.all20.B.Yp.csv")

##for H
rich.diff1=rich.diff.60m(10,H,dist,max.dist)  
dist.all1.Yp=cbind(as.data.frame(allsp),t(rich.diff1$Y),t(rich.diff1$p))
write.csv(dist.all1.Yp,"d:/rich/dist.all1.H.Yp.csv")

rich.diff5=rich.diff.60m(50,H,dist,max.dist)  
dist.all5.Yp=cbind(as.data.frame(allsp),t(rich.diff5$Y),t(rich.diff5$p))
write.csv(dist.all5.Yp,"d:/rich/dist.all5.H.Yp.csv")

rich.diff10=rich.diff.60m(100,H,dist,max.dist) 
dist.all10.Yp=cbind(as.data.frame(allsp),t(rich.diff10$Y),t(rich.diff10$p))
write.csv(dist.all10.Yp,"d:/rich/dist.all10.H.Yp.csv")

rich.diff20=rich.diff.60m(200,H,dist,max.dist) 
dist.all20.Yp=cbind(as.data.frame(allsp),t(rich.diff20$Y),t(rich.diff20$p))
write.csv(dist.all20.Yp,"d:/rich/dist.all20.H.Yp.csv")

#end of all codes.

