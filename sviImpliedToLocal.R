library(stinepack)
source("SVI_Fitting/svi.R")

TotalToLocal<-function(k,w,dwt,dwk,dwk2){
  tmp1<-(1-0.5*k/w*dwk)^2
  tmp2<--0.25*(0.25+1/w)*dwk^2
  tmp3<-0.5*dwk2
  return(dwt/(tmp1+tmp2+tmp3))
}

### get local variance grid from implied vol surface
sviImpliedToLocal<-function(sviMatrix,k_min,k_max,dk,t_min,t_max,dt){
  ### set total variance and local variance grid
  nk<-floor((k_max-k_min)/dk)+1
  nt<-floor((t_max-t_min)/dt)+1
  ks<-seq(k_min,k_max,dk)
  ts<-seq(t_min,t_max,dt)
  
  localVars<-data.frame(matrix(nrow=nt,ncol=nk))
  colnames(localVars)<-ks
  rownames(localVars)<-ts
  
  ks1<-c(k_min-dk,k_min,ks+dk)
  ts1<-c(t_min-dt,t_min,ts+dt)
  totalVars<-data.frame(matrix(nrow=nt+2,ncol=nk+2))
  colnames(totalVars)<-ks1
  rownames(totalVars)<-ts1
  
  ### get total variance grid
  ### stineman spline for each log-strike
  for(k in ks1) {
    x<-sviMatrix[,"texp"] ### input texps
    y<-numeric(length(x)) ### input total variance
    for(i in 1:length(x))
      y[i]<-svi(sviMatrix[i,2:6],k)
    x<-c(0,x)
    y<-c(0,y)
    totalVars[,toString(k)]<-stinterp(x,y,ts1)$y
  }
  
  ### get local variance grid from total variance grid
  for(i in 1:nt){
    for(j in 1:nk){
      k<-ks[j]
      w<-totalVars[i+1,j+1]
      dwt<-(totalVars[i+2,j+1]-totalVars[i,j+1])/2/dt
      dwk<-(totalVars[i+1,j+2]-totalVars[i+1,j])/2/dk
      dwk2<-(totalVars[i+1,j+2]+totalVars[i+1,j]-2*w)/dk/dk
      localVars[i,j]<-TotalToLocal(k,w,dwt,dwk,dwk2)
    }
  }
  return(localVars)
}

### get local variance at (k,t) from a local variance grid
### linear interpolation
getLocalVar0<-function(LocalGrid,k_min,dk,t_min,dt,k,t){
  nk<-round((k-k_min)/dk)+1
  nt<-round((t-t_min)/dt)+1
  return(LocalGrid[nt,nk])
}

load("sviMatrix.RData")
k_min=-0.295
k_max=0.295
dk=0.001
t_min=0.01
t_max=1.4
dt=0.01
localVar<-sviImpliedToLocal(sviMatrix,k_min,k_max,dk,t_min,t_max,dt)
#localVol<-sqrt(localVar)
time<-seq(t_min,t_max,by=dt)
strikes<-seq(k_min,k_max,by=dk)
#persp(time, strikes, as.matrix(localVar), col="green", phi=30, theta=150,ylab="Log-strike k",xlab="Time to Expiration",zlab="Local Variance", main="Local Variance Surface")
