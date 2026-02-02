#__________________________________________________________________________________________________________________
# Figure 2  ####
# latitudinal trend (all species >1 individual regardless of significant level)
#__________________________________________________________________________________________________________________
rm(list = ls())
source("myplot.r")
path = 'Sites/'
site.summary = read.csv('SiteSummary.csv',header = T)
area=site.summary$area
sname=dir(path)
j = 1           # radius unit (j=1,2,...,8 --> 2,4,...,16 m)
p.tr = 1        # set this to 0.05 for excluding non-significant species; 1 for no significant level
abund.tr = 1    # threshold excluding rare species 

R.pN = R.pS = SE.pN = SE.pS = numeric(17)
for (i in 1:17){
  cat(i,"\r")
  fname = paste0(path,sname[i],'/null.model/RNN.pN.csv')
  dat.N = read.csv(fname, fileEncoding = "GBK")
  
  fname = paste0(path,sname[i],'/null.model/RNS.pS.csv')
  dat.S = read.csv(fname, fileEncoding = "GBK")
  
  abund = dat.N[,27]
  
  use =  !is.na(dat.N[,(2+j)]) & abund>=abund.tr & dat.N[,(10+j)]<=p.tr
  n = sum(use)
  R.pN[i] = sum(dat.N[use,(2+j)]>1)/n
  SE.pN[i] = sqrt(R.pN[i]*(1-R.pN[i])/n)
  
  
  use =  !is.na(dat.N[,(2+j)]) & abund>=abund.tr & dat.S[,(10+j)]<=p.tr
  n = sum(use)
  R.pS[i] = sum(dat.S[use,(2+j)]>1)/n
  SE.pS[i] = sqrt(R.pS[i]*(1-R.pS[i])/n)
  
}

par(mfrow = c(1, 2))
myplot(x = site.summary$adj.lat.abs, R.pN, SE.pN, R.pS, SE.pS,
       xlabel = "Absolute adjusted latitude",
       ylabel = "positive neighborhood",
       c(0,70),c(0,0.7))

#__________________________________________________________________________________________________________________
# Figure 3        ####
# legumes and non-AM exclusion analysis
#__________________________________________________________________________________________________________________
rm(list = ls())
site.summary = read.csv('SiteSummary.csv',header = T)
path = 'Sites/'
sname=dir(path)
library(BSDA)

R.pN0 = R.pS0 = matrix(nrow=17,ncol=8)
R.pN1 = R.pS1 = matrix(nrow=17,ncol=8)
R.pN2 = R.pS2 = matrix(nrow=17,ncol=8)
for (i in 1:17){
  cat(i,"\r")
  fname  =paste0(path,sname[i],'/null.model/RNN.pN.csv')
  dat.N = read.csv(fname, fileEncoding = "GBK")
  
  fname  =paste0(path,sname[i],'/null.model/RNS.pS.csv')
  dat.S = read.csv(fname, fileEncoding = "GBK")
  
  abund = dat.N[,27]
  
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
  
  use0 = !is.na(dat.N[,3]) & abund>9
  use1 = !is.na(dat.N[,3]) & abund>9 & !index.LG
  use2 = !is.na(dat.N[,3]) & abund>9 & !index.FG
  
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

# bootstrappingto estimate the statistical significance of the change in R2
# observed between the full dataset and the subset following exclusions
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
  
  h1[j] = z.test(r1,r2,sigma.x = sd(r1),sigma.y = sd(r2))$p.value
  h2[j] = z.test(r1,r3,sigma.x = sd(r1),sigma.y = sd(r3))$p.value
}


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

# bootstrappingto estimate the statistical significance of the change in R2
# observed between the full dataset and the subset following exclusions
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
  
  h1[j] = z.test(r1,r2,sigma.x = sd(r1),sigma.y = sd(r2))$p.value
  h2[j] = z.test(r1,r3,sigma.x = sd(r1),sigma.y = sd(r3))$p.value
}


#__________________________________________________________________________________________________________________
# Figure 4 and Extended Data Figure 4    #####
# RNN and RNS as function of abundance for each site
#__________________________________________________________________________________________________________________
rm(list = ls())
site.summary = read.csv('SiteSummary.csv',header = T)
path = 'Sites/'
sname=dir(path)
site.name=c("Barro Colorado Island (BCI)","Baishanzu","Chebaling", "Dinghushan", "Gutianshan", "Heishiding", "Jianfengling", "Luquillo", 
         "Nanling", "Palanan", "Pasoh", "Puer", "Rabi", "TPK", "Utah", "Wanang","Yasuni")

plt=c(1:17)
par(mfrow = c(4, 2), mar=c(4.1,4,2,1))
for (i in 1:17){
  j = plt[i]
  fname  =paste0(path,sname[j],'/null.model/RNN.pN.csv')
  dat.N = read.csv(fname, fileEncoding = "GBK")
  
  fname  =paste0(path,sname[j],'/null.model/RNS.pS.csv')
  dat.S = read.csv(fname, fileEncoding = "GBK")
  
  abund = dat.N[,27]
  
  
  if (i %in% c(4,8,12,16,17) ){
    plot(dat.N[,3],abund, log='y', xlab = 'Relative neighborhood abundance', ylab = 'Abundance', xlim = c(0,2))
  }else{
    plot(dat.N[,3],abund, log='y', ylab = 'Abundance', xlab = '', xlim = c(0,2))
  }
  
  use = dat.N[,11]<0.05
  points(dat.N[use,3],abund[use],col='red')
  abline(v=1)
  title(site.name[j])
  
  if (i %in% c(4,8,12,16,17) ){
    plot(dat.S[,3],abund, log='y', xlab = 'Relative neighborhood richness', ylab='', xlim = c(0,2))
  }else{
    plot(dat.S[,3],abund, log='y', xlab = '', ylab='', xlim = c(0,2))
  }
  
  use = dat.S[,11]<0.05
  points(dat.S[use,3],abund[use],col='red')
  abline(v=1)
  title(site.name[j])
  
}

#__________________________________________________________________________________________________________________
# Figure 5     #####
# trend with mean annual temperature (MAT)
#__________________________________________________________________________________________________________________
rm(list = ls())
source("myplot.r")
path = 'Sites/'
site.summary = read.csv('SiteSummary.csv',header = T)
area=site.summary$area
sname=dir(path)
j = 1           # radius unit (j=1,2,...,8 --> 2,4,...,16 m)
R.pN = SE.pN = numeric(17)
R.pS = SE.pS = numeric(17)

for (i in 1:17){
  cat(i,"\r")
  fname = paste0(path,sname[i],'/null.model/RNN.pN.csv')
  dat.N = read.csv(fname, fileEncoding = "GBK")
  
  fname  =paste0(path,sname[i],'/null.model/RNS.pS.csv')
  dat.S = read.csv(fname, fileEncoding = "GBK")
  
  abund = dat.N[,27]
  
  use =  !is.na(dat.N[,3]) & abund>1
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



#__________________________________________________________________________________________________________________
# Extended Data Figure 1 #####
# latitudinal trend (excluding rare species, abundant<50)
#__________________________________________________________________________________________________________________
rm(list = ls())
source("myplot.r")
path = 'Sites/'
site.summary = read.csv('SiteSummary.csv',header = T)
area=site.summary$area
sname=dir(path)
j = 1           # radius unit (j=1,2,...,8 --> 2,4,...,16 m)
p.tr = 1        # set this to 0.05 for excluding non-significant species; 1 for no significant level
abund.tr = 50    # set this to 50 for excluding rare species 

R.pN = R.pS = SE.pN = SE.pS = numeric(17)
for (i in 1:17){
  cat(i,"\r")
  fname = paste0(path,sname[i],'/null.model/RNN.pN.csv')
  dat.N = read.csv(fname, fileEncoding = "GBK")
  
  fname = paste0(path,sname[i],'/null.model/RNS.pS.csv')
  dat.S = read.csv(fname, fileEncoding = "GBK")
  
  abund = dat.N[,27]
  
  use =  !is.na(dat.N[,(2+j)]) & abund>=abund.tr & dat.N[,(10+j)]<=p.tr
  n = sum(use)
  R.pN[i] = sum(dat.N[use,(2+j)]>1)/n
  SE.pN[i] = sqrt(R.pN[i]*(1-R.pN[i])/n)
  
  
  use =  !is.na(dat.N[,(2+j)]) & abund>=abund.tr & dat.S[,(10+j)]<=p.tr
  n = sum(use)
  R.pS[i] = sum(dat.S[use,(2+j)]>1)/n
  SE.pS[i] = sqrt(R.pS[i]*(1-R.pS[i])/n)
  
}

par(mfrow = c(1, 2))
myplot(x = site.summary$adj.lat.abs, R.pN, SE.pN, R.pS, SE.pS,
       xlabel = "Absolute adjusted latitude",
       ylabel = "positive neighborhood",
       c(0,70),c(0,0.7))

#__________________________________________________________________________________________________________________
# Extended Data Fig.2  #####
# latitudinal trend (excluding non-significant species according to null.model)
#__________________________________________________________________________________________________________________
rm(list = ls())
source("myplot.r")
path = 'Sites/'
site.summary = read.csv('SiteSummary.csv',header = T)
area=site.summary$area
sname=dir(path)
j = 1           # radius unit (j=1,2,...,8 --> 2,4,...,16 m)
p.tr = 0.05     # set this to 0.05 for excluding non-significant species; 1 for no significant level
abund.tr = 1    # set this to 50 for excluding rare species 

R.pN = R.pS = SE.pN = SE.pS = numeric(17)
for (i in 1:17){
  cat(i,"\r")
  fname = paste0(path,sname[i],'/null.model/RNN.pN.csv')
  dat.N = read.csv(fname, fileEncoding = "GBK")
  
  fname = paste0(path,sname[i],'/null.model/RNS.pS.csv')
  dat.S = read.csv(fname, fileEncoding = "GBK")
  
  abund = dat.N[,27]
  
  use =  !is.na(dat.N[,(2+j)]) & abund>=abund.tr & dat.N[,(10+j)]<=p.tr
  n = sum(use)
  R.pN[i] = sum(dat.N[use,(2+j)]>1)/n
  SE.pN[i] = sqrt(R.pN[i]*(1-R.pN[i])/n)
  
  
  use =  !is.na(dat.N[,(2+j)]) & abund>=abund.tr & dat.S[,(10+j)]<=p.tr
  n = sum(use)
  R.pS[i] = sum(dat.S[use,(2+j)]>1)/n
  SE.pS[i] = sqrt(R.pS[i]*(1-R.pS[i])/n)
  
}

par(mfrow = c(1, 2))
myplot(x = site.summary$adj.lat.abs, R.pN, SE.pN, R.pS, SE.pS,
       xlabel = "Absolute adjusted latitude",
       ylabel = "positive neighborhood",
       c(0,70),c(0,0.8))

#__________________________________________________________________________________________________________________
# Extended Data Fig.3      #########
# # latitudinal trend (excluding large trees)
#__________________________________________________________________________________________________________________
rm(list = ls())
source("myplot.r")
path = 'Sites/'
site.summary = read.csv('SiteSummary.csv',header = T)
area=site.summary$area
sname=dir(path)
j = 1           # radius unit (j=1,2,...,8 --> 2,4,...,16 m)
p.tr = 1        # set this to 0.05 for excluding non-significant species; 1 for no significant level
abund.tr = 50    # set this to 50 for excluding rare species 

R95.pN = R95.pS = SE95.pN = SE95.pS = numeric(17)
R90.pN = R90.pS = SE90.pN = SE90.pS = numeric(17)
for (i in 1:17){
  cat(i,"\r")
  fname = paste0(path,sname[i],'/null.model/RNN95.pN.csv')
  dat.N = read.csv(fname, fileEncoding = "GBK")
  
  fname = paste0(path,sname[i],'/null.model/RNS95.pS.csv')
  dat.S = read.csv(fname, fileEncoding = "GBK")
  
  abund = dat.N[,11]
  
  use =  !is.na(dat.N[,(2+j)]) & abund>=abund.tr
  n = sum(use)
  R95.pN[i] = sum(dat.N[use,(2+j)]>1)/n
  SE95.pN[i] = sqrt(R95.pN[i]*(1-R95.pN[i])/n)
  
  
  use =  !is.na(dat.S[,(2+j)]) & abund>=abund.tr
  n = sum(use)
  R95.pS[i] = sum(dat.S[use,(2+j)]>1)/n
  SE95.pS[i] = sqrt(R95.pS[i]*(1-R95.pS[i])/n)
  
  #  90 
  fname = paste0(path,sname[i],'/null.model/RNN90.pN.csv')
  dat.N = read.csv(fname, fileEncoding = "GBK")
  
  fname = paste0(path,sname[i],'/null.model/RNS90.pS.csv')
  dat.S = read.csv(fname, fileEncoding = "GBK")
  
  abund = dat.N[,11]
  
  use =  !is.na(dat.N[,(2+j)]) & abund>=abund.tr
  n = sum(use)
  R90.pN[i] = sum(dat.N[use,(2+j)]>1)/n
  SE90.pN[i] = sqrt(R90.pN[i]*(1-R90.pN[i])/n)
  
  
  use =  !is.na(dat.S[,(2+j)]) & abund>=abund.tr
  n = sum(use)
  R90.pS[i] = sum(dat.S[use,(2+j)]>1)/n
  SE90.pS[i] = sqrt(R90.pS[i]*(1-R90.pS[i])/n)
  
}

par(mfrow = c(2, 2))
myplot(x = site.summary$adj.lat.abs, R95.pN, SE95.pN, R95.pS, SE95.pS,
       xlabel = "Absolute adjusted latitude",
       ylabel = "positive neighborhood",
       c(0,70),c(0,0.7))

myplot(x = site.summary$adj.lat.abs, R90.pN, SE90.pN, R90.pS, SE90.pS,
       xlabel = "Absolute adjusted latitude",
       ylabel = "positive neighborhood",
       c(0,70),c(0,1))


#__________________________________________________________________________________________________________________
# Extended Data Figure 6    #####
# latitudinal trend with absolute proportions
#__________________________________________________________________________________________________________________
rm(list = ls())
source("myplot.r")
site.summary = read.csv('SiteSummary.csv',header = T)
path = 'Sites/'
sname=dir(path)
j = 1           # radius unit (j=1,2,...,8 --> 2,4,...,16 m)

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

N = numeric(17)
for (i in 1:17){
  cat(i,"\r")
  fname  =paste0(path,sname[i],'/null.model/RNN.pN.csv')
  dat.N = read.csv(fname, fileEncoding = "GBK")
  
  fname  =paste0(path,sname[i],'/null.model/RNS.pS.csv')
  dat.S = read.csv(fname, fileEncoding = "GBK")
  
  abund = dat.N[,27]
  
  
  use =  !is.na(dat.N[,3]) & abund>1
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
  
  
  
  N[i] = mean(abund[use]) #mean site abundance
  
  
  ## standardized
  nd = data.frame(x  = 3)
  x = log10(abund[use])
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



#__________________________________________________________________________________________________________________
# Extended Data Figure 7 ####
# z-score
#__________________________________________________________________________________________________________________
rm(list = ls())
site.name=c("Barro Colorado Island (BCI)","Baishanzu","Chebaling", "Dinghushan", "Gutianshan", "Heishiding", "Jianfengling", "Luquillo", 
            "Nanling", "Palanan", "Pasoh", "Puer", "Rabi", "TPK", "Utah", "Wanang","Yasuni")
path = 'Sites/'
sname=dir(path)
j = 1           # radius unit (j=1,2,...,8 --> 2,4,...,16 m)
par(mfrow = c(4, 2), mar=c(4.1,4,2,1))
for (i in 1:17){
  
  fname  =paste0(path,sname[i],'/null.model/RNN.pN.csv')
  dat.N = read.csv(fname, fileEncoding = "GBK")
  
  fname  =paste0(path,sname[i],'/null.model/RNS.pS.csv')
  dat.S = read.csv(fname, fileEncoding = "GBK")
  
  abund = dat.N[,27]
  
  
  if (i %in% c(4,8,12,16,17)){
    plot(dat.N[,19],abund, log='y', xlab = 'z-score (abundance)', ylab = 'Abundance', xlim = c(-20,20))
  }else{
    plot(dat.N[,19],abund, log='y', ylab = 'abundance', xlab = '', xlim = c(-20,20))
  }
  
  use = dat.N[,11]<0.05
  points(dat.N[use,19],abund[use],col='red')
  abline(v =+1.96, col = "black", lty = 2, lwd = 1)
  abline(v =-1.96, col = "black", lty = 2, lwd = 1)
  title(site.name[i])
  
  if (i %in% c(4,8,12,16,17)){
    plot(dat.S[,19],abund, log='y', xlab = 'z-score (richness)', ylab='', xlim = c(-20,20))
  }else{
    plot(dat.S[,19],abund, log='y', xlab = '', ylab='', xlim = c(-20,20))
  }
  
  use = dat.S[,11]<0.05
  points(dat.S[use,19],abund[use],col='red')
  abline(v =+1.96, col = "black", lty = 2, lwd = 1)
  abline(v =-1.96, col = "black", lty = 2, lwd = 1)
  title(site.name[i])
  
}

