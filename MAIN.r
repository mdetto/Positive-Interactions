site.summary = read.csv('SiteSummary17.csv',header = T)
area=site.summary$area
path = 'Sites/'
sname=dir(path)
library(BSDA)


###########     Fig. 2 (or extended Data Fig.2)  ###########################
j = 1  # radius unit (j=1:8)
R.pN = R.pS = numeric(17)
SE.pN = SE.pS = numeric(17)

for (i in 1:17){
  cat(i,"\r")
  fname = paste0(path,sname[i],'/null.model.matteo/RNN.pN.csv')
  dat.N = read.csv(fname, fileEncoding = "GBK")
  
  fname = paste0(path,sname[i],'/null.model.matteo/RNS.pS.csv')
  dat.S = read.csv(fname, fileEncoding = "GBK")
  
  fname = paste0(path,sname[i],'/all.individuals/No.stems1.csv')
  abund = read.csv(fname, fileEncoding = "GBK")
  
  use =  !is.na(dat.N[,3]) & abund[,2]>=1
  
  n = sum(use)
  R.pN[i] = sum(dat.N[use,(2+j)]>1)/n
  SE.pN[i] = sqrt(R.pN[i]*(1-R.pN[i])/n)
  
  R.pS[i] = sum(dat.S[use,(2+j)]>1)/n
  SE.pS[i] = sqrt(R.pS[i]*(1-R.pS[i])/n)
  
  
}


par(mfrow = c(1, 2))
myplot(x = site.summary$adj.lat.abs, R.pN, SE.pN, R.pS, SE.pS,
       xlabel = "Absolute adjusted latitude",
       ylabel = "positive neighborhood",
       c(0,70),c(0,0.7))


##############Fig. 3:legumes and non-AM exclusion analysis###############
site.summary = read.csv('SiteSummary17.csv',header = T)
path = 'Sites/'
sname=dir(path)
library(BSDA)

## 60 m, DBH class NULL model test
R.pN0 = R.pS0 = matrix(nrow=17,ncol=8)
R.pN1 = R.pS1 = matrix(nrow=17,ncol=8)
R.pN2 = R.pS2 = matrix(nrow=17,ncol=8)

for (i in 1:17){
  cat(i,"\r")
  fname  =paste0(path,sname[i],'/null.model.matteo/RNN.pN.csv')
  dat.N = read.csv(fname, fileEncoding = "GBK")
  
  fname  =paste0(path,sname[i],'/null.model.matteo/RNS.pS.csv')
  dat.S = read.csv(fname, fileEncoding = "GBK")
  
  fname  =paste0(path,sname[i],'/all.individuals/No.stems1.csv')
  abund = read.csv(fname, fileEncoding = "GBK")
  
  # legume
  fname  =paste0(path,sname[i],'/all.individuals/legumes.csv')
  legumes = read.csv(fname, fileEncoding = "GBK")  
  sp = dat.S[,2]
  index.LG = logical(length(sp))
  index.FX = logical(length(sp))
  
  for (j in 1:length(legumes[,2])){
    index.LG[sp==legumes[j,2]]=TRUE
    index.FX[sp==legumes[j,2]]= legumes[j,3]=="y"
  }
  
  # fungi
  fname  =paste0(path,sname[i],'/all.individuals/fungi.csv')
  fungi = read.csv(fname, fileEncoding = "GBK")  
  sp = dat.S[,2]
  index.FG = logical(length(sp))
  
  for (j in 1:length(fungi[,2])){
    index.FG[sp==fungi[j,2] & fungi[j,3]=="AM"]=TRUE
  }
  
  use0 = !is.na(dat.N[,3]) & abund[,2]>9
  use1 = !is.na(dat.N[,3]) & abund[,2]>9 & !index.LG
  use2 = !is.na(dat.N[,3]) & abund[,2]>9 & !index.FG
  
  for (j in 1:8){
    # all species
    R.pN0[i,j] = sum(dat.N[use0,(2+j)]>1)/sum(dat.N[use0,(2+j)])
    R.pS0[i,j] = sum(dat.S[use0,(2+j)]>1)/sum(dat.S[use0,(2+j)])
    
    # legume excluded
    R.pN1[i,j] = mean(dat.N[use1,(2+j)])
    R.pS1[i,j] = mean(dat.S[use1,(2+j)])
    
    # EcM excluded
    R.pN2[i,j] = mean(dat.N[use2,(2+j)])
    R.pS2[i,j] = mean(dat.S[use2,(2+j)])
    
  }
  
}


par(mfrow = c(1, 2))

#abundance
r1=r2=r3=numeric(8)
for (i in 1:8){
  
  r1[i]=summary(lm(R.pN0[,i]~site.summary$adj.lat.abs))$coefficients[2]
  r2[i]=summary(lm(R.pN1[,i]~site.summary$adj.lat.abs))$coefficients[2]
  r3[i]=summary(lm(R.pN2[,i]~site.summary$adj.lat.abs))$coefficients[2]
  
}

plot(c(1:8)*2,r1,type='o',col = 'black',xlab = 'Distance to the focal individual (m)',ylab = 'Slope between relative neighborhood abundance and latitude',ylim=c(-0.008,0.002),pch=19)
points(c(1:8)*2,r2,col='red',pch=19)
lines(c(1:8)*2,r2,col='red')
points(c(1:8)*2,r3,col='blue',pch=19)
lines(c(1:8)*2,r3,col='blue')
legend(x='bottomright',legend=c('All species','Leguminous species excluded',"Non-AM species excluded"), fill = c("black","red", "blue") )

## bootstrapping to get the p values, for the before and after exclusion comparison
r1=r2=numeric(999)
h1 = h2 = numeric(8)
for (j in 1:8){
  cat(j,'')
  for (i in 1:999){
    
    use = sample(17, replace = TRUE);
    r1[i]=summary(lm(R.pN0[use,j]~site.summary$adj.lat.abs[use]))$coefficients[2]  
    r2[i]=summary(lm(R.pN1[use,j]~site.summary$adj.lat.abs[use]))$coefficients[2]
    r3[i]=summary(lm(R.pN2[use,j]~site.summary$adj.lat.abs[use]))$coefficients[2]
    
  }
  
  # h[j] = t.test(r1,r2)$p.value
  h1[j] = z.test(r1,r2,sigma.x = sd(r1),sigma.y = sd(r2))$p.value
  h2[j] = z.test(r1,r3,sigma.x = sd(r1),sigma.y = sd(r3))$p.value
}
h1
h2
as.numeric(format(round(h1, 4)))
format(round(h2, 4),nsmall = 4)

# richness
lR.pH0=log(R.pS0)
lR.pH2=log(R.pS2)
site.summary$ladj.lat.abs=log(site.summary$adj.lat.abs)
r1=r2=r3=numeric(8)
for (i in 1:8){
  
  r1[i]=summary(lm(R.pS0[,i]~site.summary$adj.lat.abs))$coefficients[2]
  r2[i]=summary(lm(R.pS1[,i]~site.summary$adj.lat.abs))$coefficients[2]
  r3[i]=summary(lm(R.pS2[,i]~site.summary$adj.lat.abs))$coefficients[2]
  
}

plot(c(1:8)*2,r1,type='o',col = 'black',xlab = 'Distance to the focal individual (m)',ylab = 'Slope between relative neighborhood richness and latitude',ylim=c(-0.008,0.002),pch=19)
points(c(1:8)*2,r2,col='red',pch=19)
lines(c(1:8)*2,r2,col='red')
points(c(1:8)*2,r3,col='blue',pch=19)
lines(c(1:8)*2,r3,col='blue')
legend(x='bottomright',legend=c('All species','Leguminous species excluded',"Non-AM species excluded"), fill = c("black","red", "blue") )

## bootstrapping to get the p values, for the before and after exclusion comparison
r1=r2=numeric(999)
h1 = h2 = numeric(8)
for (j in 1:8){
  cat(j,'')
  for (i in 1:999){
    
    use = sample(17, replace = TRUE);
    r1[i]=summary(lm(R.pS0[use,j]~site.summary$adj.lat.abs[use]))$coefficients[2]  
    r2[i]=summary(lm(R.pS1[use,j]~site.summary$adj.lat.abs[use]))$coefficients[2]
    r3[i]=summary(lm(R.pS2[use,j]~site.summary$adj.lat.abs[use]))$coefficients[2]
    
  }
  
  # h[j] = t.test(r1,r2)$p.value
  h1[j] = z.test(r1,r2,sigma.x = sd(r1),sigma.y = sd(r2))$p.value
  h2[j] = z.test(r1,r3,sigma.x = sd(r1),sigma.y = sd(r3))$p.value
}
h1
h2
as.numeric(format(round(h1, 4)))
format(round(h2, 4),nsmall = 4)



###########  Fig. 5   ###########################
j = 1
R.pN = R.pS = numeric(17)
SE.pN = SE.pS = numeric(17)

for (i in 1:17){
  cat(i,"\r")
  fname = paste0(path,sname[i],'/null.model.matteo/RNN.pN.csv')
  dat.N = read.csv(fname, fileEncoding = "GBK")
  
  fname  =paste0(path,sname[i],'/null.model.matteo/RNS.pS.csv')
  dat.S = read.csv(fname, fileEncoding = "GBK")
  
  fname  =paste0(path,sname[i],'/all.individuals/No.stems1.csv')
  abund = read.csv(fname, fileEncoding = "GBK")
  
  use =  !is.na(dat.N[,3]) & abund[,2]>1
  n = sum(use)
  
  R.pN[i] = sum(dat.N[use,(2+j)]>1)/n
  SE.pN[i] = sqrt(R.pN[i]*(1-R.pN[i])/n)
  
  R.pS[i] = sum(dat.S[use,(2+j)]>1)/n
  SE.pS[i] = sqrt(R.pS[i]*(1-R.pS[i])/n)
  
}


par(mfrow = c(1, 2))
myplot(site.summary$annu.temp, R.pN, SE.pN, R.pS, SE.pS,
       xlabel = "MAT (C)",
       ylabel = "positive neighborhood",
       c(0,30),c(0,0.7))



###########  Extended Data Figure 1 ###########################
j = 1
R.pN = R.pS = numeric(17)
SE.pN = SE.pS = N= numeric(17)

for (i in 1:17){
  cat(i,"\r")
  fname  =paste0(path,sname[i],'/null.model.matteo/RNN.pN.csv')
  dat.N = read.csv(fname, fileEncoding = "GBK")
  
  fname  =paste0(path,sname[i],'/null.model.matteo/RNS.pS.csv')
  dat.S = read.csv(fname, fileEncoding = "GBK")
  
  fname  =paste0(path,sname[i],'/all.individuals/No.stems1.csv')
  abund = read.csv(fname, fileEncoding = "GBK")
  
  use =  !is.na(dat.N[,3]) & abund[,2]>=1
  
  n = sum(dat.N[use,(10+j)]<0.05)
  R.pN[i] = sum(dat.N[use,(2+j)]>1 & dat.N[use,(10+j)]<0.05)/n
  SE.pN[i] = sqrt(R.pN[i]*(1-R.pN[i])/n)
  
  n = sum(dat.S[use,(10+j)]<0.05)
  R.pS[i] = sum(dat.S[use,(2+j)]>1 & dat.S[use,(10+j)]<0.05)/n
  SE.pS[i] = sqrt(R.pS[i]*(1-R.pS[i])/n)
  
  N[i] = mean(as.numeric(abund[use,2]))
}


par(mfrow = c(1, 2))
j = 1

myplot(x = site.summary$adj.lat.abs, R.pN, SE.pN, R.pS, SE.pS,
       xlabel = "Absolute adjusted latitude",
       ylabel = "positive neighborhood",
       c(0,70),c(0,0.7))


###########      Extended data Fig. 2    ###########################
par(mfrow = c(2, 2))

plot(site.summary$adj.lat.abs,R.pN1[,1],xlab = "Absolute adjusted latitude", ylab = c("Proportion of species with","positive neighborhood abundance"),
     xlim = c(0,70), ylim = c(0.07,0.53), pch=19)
lm1 = lm(R.pN1[,1]~site.summary$adj.lat.abs)
abline(lm1)
z=summary(lm1)
r2 = z$r.square
p = z$coefficients[8]

if (p<1e-3){
  title(paste0('a) RL: R2 = ',as.character(round(r2,2)),', p<1e-3'))
} else {
  title(paste0('a) RL: R2 = ',as.character(round(r2,2)),', p  = ',as.character(round(p,4))))
}


plot(site.summary$adj.lat.abs,R.pN2[,1],xlab = "Absolute adjusted latitude", ylab = c("Proportion of species with","positive neighborhood abundance"),
     xlim = c(0,70), ylim = c(0.7,1), pch=19)
lm1 = lm(R.pN2[,1]~site.summary$adj.lat.abs)
abline(lm1)
z=summary(lm1)
r2 = z$r.square
p = z$coefficients[8]
if (p<1e-3){
  title(paste0('b) CSR: R2 = ',as.character(round(r2,2)),', p<1e-3'))
} else {
  title(paste0('b) CSR: R2 = ',as.character(round(r2,2)),', p  = ',as.character(round(p,4))))
}

## removing Utah plot
abline(lm(R.pN2[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17),1]~site.summary$adj.lat.abs[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17)]),col='red')
summary(lm(R.pN2[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17),1]~site.summary$adj.lat.abs[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17)]))


plot(site.summary$adj.lat.abs,R.pS1[,1],xlab = "Absolute adjusted latitude", ylab = c("Proportion of species with","positive neighborhood richness"),
     xlim = c(0,70), ylim = c(0.,0.2), pch=19)
lm1 = lm(R.pS1[,1]~site.summary$adj.lat.abs)
abline(lm1)
z=summary(lm1)
r2 = z$r.square
p = z$coefficients[8]
if (p<1e-3){
  title(paste0('c) RL: R2 = ',as.character(round(r2,2)),', p<1e-3'))
} else {
  title(paste0('c) RL: R2 = ',as.character(round(r2,2)),', p  = ',as.character(round(p,4))))
}


plot(site.summary$adj.lat.abs,R.pS2[,1],xlab = "Absolute adjusted latitude", ylab = c("Proportion of species with","positive neighborhood richness"),
     xlim = c(0,70), ylim = c(0,1), pch=19)
lm1 = lm(R.pS2[,1]~site.summary$adj.lat.abs)
abline(lm1)
z=summary(lm1)
r2 = z$r.square
p = z$coefficients[8]
if (p<1e-3){
  title(paste0('d) CSR: R2 = ',as.character(round(r2,2)),', p<1e-3'))
} else {
  title(paste0('d) CSR: R2 = ',as.character(round(r2,2)),', p  = ',as.character(round(p,4))))
}


## removing Utah plot
abline(lm(R.pS2[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17),1]~site.summary$adj.lat.abs[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17)]),col='red')
summary(lm(R.pS2[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17),1]~site.summary$adj.lat.abs[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17)]))



sname1=c("Barro Colorado Island (BCI)","Baishanzu","Chebaling", "Dinghushan", "Gutianshan", "Heishiding", "Jianfengling", "Luquillo", 
         "Nanling", "Palanan", "Pasoh", "Puer", "Rabi", "TPK", "Utah", "Wanang","Yasuni")
#########  Figure 5  ###################
plt=c(11,1,7,5)
par(mfrow = c(4, 2), mar=c(4.1,4,2,1))
for (i in 1:4){
  j = plt[i]
  fname  =paste0(path,sname[j],'/null.model.matteo/RNN.pN.csv')
  dat.N = read.csv(fname, fileEncoding = "GBK")
  
  fname  =paste0(path,sname[j],'/null.model.matteo/RNS.pS.csv')
  dat.S = read.csv(fname, fileEncoding = "GBK")
  
  fname  =paste0(path,sname[j],'/all.individuals/No.stems1.csv')
  abund = read.csv(fname, fileEncoding = "GBK")
  
  
  if (i == 4){
    plot(dat.N[,3],abund[,2], log='y', xlab = 'Relative neighborhood abundance', ylab = 'Abundance', xlim = c(0,2))
  }else{
    plot(dat.N[,3],abund[,2], log='y', ylab = 'Abundance', xlab = '', xlim = c(0,2))
  }
  
  use =     dat.N[,11]<0.05
  points(dat.N[use,3],abund[use,2],col='red')
  abline(v=1)
  title(sname1[j])
  
  if (i == 4){
    plot(dat.S[,3],abund[,2], log='y', xlab = 'Relative neighborhood richness', ylab='', xlim = c(0,2))
  }else{
    plot(dat.S[,3],abund[,2], log='y', xlab = '', ylab='', xlim = c(0,2))
  }
  
  use =     dat.S[,11]<0.05
  points(dat.S[use,3],abund[use,2],col='red')
  abline(v=1)
  title(sname1[j])
  
}

###########  Extended Figure 4  #################
plt=c(1:17)
par(mfrow = c(4, 2), mar=c(4.1,4,2,1))
for (i in 1:17){
  j = plt[i]
  fname  =paste0(path,sname[j],'/null.model.matteo/RNN.pN.csv')
  dat.N = read.csv(fname, fileEncoding = "GBK")
  
  fname  =paste0(path,sname[j],'/null.model.matteo/RNS.pS.csv')
  dat.S = read.csv(fname, fileEncoding = "GBK")
  
  fname  =paste0(path,sname[j],'/all.individuals/No.stems1.csv')
  abund = read.csv(fname, fileEncoding = "GBK")
  
  
  if (i %in% c(4,8,12,16,17) ){
    plot(dat.N[,3],abund[,2], log='y', xlab = 'Relative neighborhood abundance', ylab = 'Abundance', xlim = c(0,2))
  }else{
    plot(dat.N[,3],abund[,2], log='y', ylab = 'Abundance', xlab = '', xlim = c(0,2))
  }
  
  use =     dat.N[,11]<0.05
  points(dat.N[use,3],abund[use,2],col='red')
  abline(v=1)
  title(sname1[j])
  
  if (i %in% c(4,8,12,16,17) ){
    plot(dat.S[,3],abund[,2], log='y', xlab = 'Relative neighborhood richness', ylab='', xlim = c(0,2))
  }else{
    plot(dat.S[,3],abund[,2], log='y', xlab = '', ylab='', xlim = c(0,2))
  }
  
  use =     dat.S[,11]<0.05
  points(dat.S[use,3],abund[use,2],col='red')
  abline(v=1)
  title(sname1[j])
  
}




#----------------------------------------------------------------------------------
for (iii in 1:17){
  
  fname  =paste0('Sites/',sname[iii],'/plotname.all.Rdata')
  dat = load(fname)
  
  res$cvb=cvb
  res$Lx=Lx
  res$Ly=Ly
  res$n0=length(res$N[1,])
  fname  =paste0('Sites/',sname[iii],'/NewTest.RData')
  saveRDS(res,fname)
}

##### for LARGE TREE EXCLUSION ANALYSIS ############
site.summary = read.csv('SiteSummary17.csv',header = T)
path = 'Sites/'
sname=dir(path)

n.class = 25

for (iii in 1:17){
  cat(iii,"\r")
  fname  =paste0(path,sname[iii],'/plotname.all.RData')
  datx = load(fname)
  
  
  q = quantile(log(dbh),probs=seq(0,1,length.out = n.class+1))
  X.S=X.H=numeric(n.class)
  for (i in 1:n.class){
    use=log(dbh)>=q[i] & log(dbh)<q[i+1]
    X.S[i] = mean(res$S[1,use])
    X.H[i] = mean(res$H[1,use])
    
  }
  
  par(mfrow = c(1, 2))
  
  plot(1:n.class, X.S, type = "o") 
  plot(1:n.class, X.H/X.S, type = "o")  
  title(sname[iii])
  
}

##########  Extended Data Figure 6 ###########################
j = 1
R.aN = R.aS = numeric(17)
R.pN = R.pS = numeric(17)
R.nN = R.nS = numeric(17)
SE.aN = SE.aS = numeric(17)
SE.pN = SE.pS = numeric(17)
SE.nN = SE.nS = numeric(17)

R.pN.std = R.pS.std = numeric(17)
R.nN.std = R.nS.std = numeric(17)
SE.pN.std = SE.pS.std = numeric(17)
SE.nN.std = SE.nS.std = numeric(17)


for (i in 1:17){
  cat(i,"\r")
  fname  =paste0(path,sname[i],'/null.model.matteo/RNN.pN.csv')
  dat.N = read.csv(fname, fileEncoding = "GBK")
  
  fname  =paste0(path,sname[i],'/null.model.matteo/RNS.pS.csv')
  dat.S = read.csv(fname, fileEncoding = "GBK")
  
  fname  =paste0(path,sname[i],'/all.individuals/No.stems1.csv')
  abund = read.csv(fname, fileEncoding = "GBK")
  
  
  use =  !is.na(dat.N[,3]) & abund[,2]>1
  n = sum(use)
  
  R.pN[i] = sum(dat.N[use,(2+j)]>1 & dat.N[use,(10+j)]<0.05)/n
  R.pS[i] = sum(dat.S[use,(2+j)]>1 & dat.S[use,(10+j)]<0.05)/n
  SE.pN[i] = sqrt(R.pN[i]*(1-R.pN[i])/n)
  SE.pS[i] = sqrt(R.pS[i]*(1-R.pS[i])/n)
  
  R.nN[i] = sum(dat.N[use,(2+j)]<1 & dat.N[use,(10+j)]<0.05)/n
  R.nS[i] = sum(dat.S[use,(2+j)]<1 & dat.S[use,(10+j)]<0.05)/n
  SE.nN[i] = sqrt(R.nN[i]*(1-R.nN[i])/n)
  SE.nS[i] = sqrt(R.nS[i]*(1-R.nS[i])/n)
  
  R.aN[i] = sum(dat.N[use,(10+j)]<0.05)/n
  R.aS[i] = sum(dat.S[use,(10+j)]<0.05)/n
  SE.aN[i] = sqrt(R.aN[i]*(1-R.aN[i])/n)
  SE.aS[i] = sqrt(R.aS[i]*(1-R.aS[i])/n)
  
  
  
  N[i] = mean(as.numeric(abund[use,2]))
  
  
  ## standardized
  nd = data.frame(x  = 3)
  x = log10(abund[use,2])
  y = dat.N[use,(2+j)]>1 & dat.N[use,(10+j)]<0.05  
  lm1 = glm(y~x,family=binomial())
  pred = predict(lm1,nd,type = "response", se.fit = TRUE)
  R.pN.std[i]  = pred$fit
  SE.pN.std[i] = pred$se.fit
  
  y = dat.N[use,(2+j)]>1 & dat.S[use,(10+j)]<0.05  
  lm1=glm(y~x,family=binomial())
  pred = predict(lm1,nd,type = "response", se.fit = TRUE);
  R.pS.std[i]  = pred$fit
  SE.pS.std[i] = pred$se.fit
  
  y = dat.N[use,(2+j)]<1 & dat.N[use,(10+j)]<0.05  
  lm1=glm(y~x,family=binomial())
  pred = predict(lm1,nd,type = "response", se.fit = TRUE);
  R.nN.std[i]  = pred$fit
  SE.nN.std[i] = pred$se.fit
  
  y = dat.N[use,(2+j)]<1 & dat.S[use,(10+j)]<0.05  
  lm1=glm(y~x,family=binomial())
  pred = predict(lm1,nd,type = "response", se.fit = TRUE);
  R.nS.std[i]  = pred$fit
  SE.nS.std[i] = pred$se.fit
  
}




par(mfrow = c(2, 2))
## positive
myplot(x = site.summary$adj.lat.abs, R.pN, SE.pN, R.pS, SE.pS,
       xlabel = "Absolute adjusted latitude",
       ylabel = "positive neighborhood",
       c(0,70),c(0,1))



## negative
myplot(x = site.summary$adj.lat.abs, R.nN, SE.nN, R.nS, SE.nS,
       xlabel = "Absolute adjusted latitude",
       ylabel = "negative neighborhood",
       c(0,70),c(0,1))

#####
par(mfrow = c(1, 2))
## all significant
myplot(x = N, R.aN, SE.aN, R.aS, SE.aS,
       xlabel = "mean species abundance",
       ylabel = "significant neighborhood",
       c(100,2000), c(0,1))


### standardized
par(mfrow = c(2, 2))
## positive
myplot(x = site.summary$adj.lat.abs, R.pN.std, SE.pN.std*0, R.pS.std, SE.pS.std*0,
       xlabel = "Absolute adjusted latitude",
       ylabel = "positive neighborhood",
       c(0,70), c(0,1))

## negative
myplot(x = site.summary$adj.lat.abs, R.nN.std, SE.nN.std*0, R.nS.std, SE.nS.std*0,
       xlabel = "Absolute adjusted latitude",
       ylabel = "negative neighborhood",
       c(0,70), c(0,1))



