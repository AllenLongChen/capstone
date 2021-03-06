setwd("/Users/mengwang/Documents/Capstone")
source("BlackScholes.R")
source("localVolMC.R")
source("localVolPDE.R")
source("localVolDupirePDE.R")
source("LRBinomial.R")

### set local vol parameters
params <- list(a=0.0012,b=0.1634,sig=0.1029,rho=-0.5555,m1=0.0439)
#params <- list(a=0.04,b=0,sig=0.1029,rho=-0.5555,m1=0.0439)

### local vol MC
#strikes <- c(0.8,1.0,1.2)
strikes<-seq(0.5,1.5,0.02)
callPriceMC<-LocalVolMC(params)(S0=1, capT=1, strikes, N=10000, m=16, evolve=evolveEuler, exactVols=NULL)
impliedVolMC<-BSImpliedVolCall(S0=1, K=strikes, capT=1, r=0, C=callPriceMC)
plot(log(strikes),impliedVolMC,ylim=c(0.1,0.4),ylab="impliedVol",xlab="log-strike")

### local vol Dupire PDE
sigmaLocalVolGeneric <- function(params) {
  a<-params$a
  b<-params$b
  sig<-params$sig
  rho<-params$rho
  m1<-params$m1
  res <- function(k,t){
    v<-a+b*(rho*(k/sqrt(t)-m1)+sqrt((k/sqrt(t)-m1)^2+sig^2*t))
    return(sqrt(v))
  }
  return(res)
}

sigmaLocalVol<-sigmaLocalVolGeneric(params)

callPricePDE<-numeric(length(strikes))
impliedVolPDE<-numeric(length(strikes))
callPriceLR<-numeric(length(strikes))
impliedVolLR<-numeric(length(strikes))
for (i in 1:length(strikes)){
  strike=strikes[i]
  #callPricePDE[i]<-callLocalVolDupirePDE(S0=1,K=strike,r=0, q=0, sigmaLocalVol, t=1, dK=0.01, dt=0.01, sdw=4, start=1, theta=1/2)
  callPricePDE[i]<-callLocalVolPDE(S0=1,K=strike,r=0, q=0, sigmaLocalVol, t=1, dS=0.01, dt=0.01, sdw=4, start=1, theta=1/2)
  impliedVolPDE[i]<-BSImpliedVolCall(S0=1, K=strike, capT=1, r=0, C=callPricePDE[i])
  callPriceLR[i]<-binPriceLR(X=strike, S0=1, r=0, sig=0.2, Texp=1, oType="C", earlyExercise=FALSE, steps=100)
  #impliedVolLR[i]<-binImpVolLR(S0=1, K=strike, capT=1, r=0, V=callPriceLR[i], oType="C", earlyExercise=FALSE, steps=10)
}

#print(callPricePDE)
#print(impliedVolPDE)
#print(callPriceLR)
#print(impliedVolLR)

#print(BSFormula(S0=1, K=strikes, capT=1, r=0, sigma=0.2))
lines(log(strikes),impliedVolPDE,col="red")

