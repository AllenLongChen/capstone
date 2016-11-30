source("eurToAmerLV.R")
source("localVolPDEAmerican.R")

load("GOOG.RData")

df[,"mid"]<-(df[,"best_bid"]+df[,"best_offer"])/2

df<-df[,c("strike_price","Texp","mid","cp_flag")]
df[,"mid"]<-df[,"mid"]/869.81 
df[,"strike_price"]<-df[,"strike_price"]/869.81

df_c<-subset(df,cp_flag=="C")
df_p<-subset(df,cp_flag=="P")

#svi_param_put<-list(a=0.03030448,b=0.12966745,sigma=0.21952271,rho=-0.54106085,m=-0.05339763)
#svi_param_call<-list(a=0.03260099,b=0.02364996,sigma=0.48849012,rho=-0.55487621,m=0.11030024)

svi_param_put<-list(a=0.029721295,b=0.081014508,sigma=0.209072884,rho=-0.567375142,m=0.054033291)
svi_param_call<-list(a=0.029721295,b=0.081014508,sigma=0.209072884,rho=-0.567375142,m=0.054033291)

rs<-c(-0.004284978,-0.001918212,-0.013127719,-0.036206437,-0.007581833,
      0.006165864,0.000614612,-0.00014151,0.001334115,0.002855994)

sigmaLocalVolGeneric <- function(params) {
  a<-params$a
  b<-params$b
  sig<-params$sig
  rho<-params$rho
  m<-params$m
  res <- function(k,t){
    v<-a+b*(rho*(k/sqrt(t)-m)+sqrt((k/sqrt(t)-m)^2+sig^2*t))
    return(sqrt(v))
  }
  return(res)
}

sigmaLocalVol_put<-sigmaLocalVolGeneric(svi_param_put)
sigmaLocalVol_call<-sigmaLocalVolGeneric(svi_param_call)

Texps<-unique(df[,"Texp"])

### put option
par(mfrow=c(2,5),mex=0.5)
for(i in 1:length(Texps)){
  print(i)
  texp<-Texps[i]
  df_p1<-subset(df_p,Texp==texp)
  params_put <- list(S0 = 1, r = rs[i], iscall = FALSE)
  ameprice_put <- getAmericanPriceFunc(params_put, SVI_LocalVol_func(svi_param_put,params_put$S0), EurImpVol_func(svi_param_put, params_put$S0))
  ap_lvol<-ameprice_put(texp,df_p1[,"strike_price"])
  ap_pde<-sapply(df_p1[,"strike_price"],function(x){LocalVolPDEAmerican(S0=1, K=x, r=rs[i], q=0,
                      sigmaLocalVol_put,t=texp, dS=0.01,dt=texp/100, sdw=4, start=1,
                      theta=1/2,oType="P")})
  title=paste("T=",round(texp*1000)/1000,sep="")
  plot(log(df_p1[,"strike_price"]),ap_lvol,type="l",col="red",xlab="Log-Strike",ylab="American Put",main=title)
  lines(log(df_p1[,"strike_price"]),ap_pde,col="blue")
  lines(log(df_p1[,"strike_price"]),df_p1[,"mid"])
  #legend("topleft",c("early boundary approx","Dupire PDE","mid price"),col=c("red","blue","black"))
}

### call option
par(mfrow=c(2,5),mex=0.5)
for(i in 1:length(Texps)){
  print(i)
  texp<-Texps[i]
  df_c1<-subset(df_c,Texp==texp)
  params_call <- list(S0 = 1, r = rs[i], iscall = TRUE)
  ameprice_call <- getAmericanPriceFunc(params_call, SVI_LocalVol_func(svi_param_call,params_call$S0), EurImpVol_func(svi_param_call, params_call$S0))
  ap_lvol<-ameprice_call(texp,df_c1[,"strike_price"])
  ap_pde<-sapply(df_c1[,"strike_price"],function(x){LocalVolPDEAmerican(S0=1, K=x, r=rs[i], q=0,
                                                                        sigmaLocalVol_call,t=texp, dS=0.01,dt=texp/100, sdw=4, start=1,
                                                                        theta=1/2,oType="C")})
  title=paste("T=",round(texp*1000)/1000,sep="")
  plot(log(df_c1[,"strike_price"]),ap_lvol,type="l",col="red",xlab="Log-Strike",ylab="American Call",main=title)
  lines(log(df_c1[,"strike_price"]),ap_pde,col="blue")
  lines(log(df_c1[,"strike_price"]),df_c1[,"mid"])
  #legend("topleft",c("early boundary approx","Dupire PDE","mid price"),col=c("red","blue","black"))
}