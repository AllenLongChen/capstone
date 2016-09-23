source("BlackScholes.R")
source("localVolMC.R")
source("localVolDupirePDE.R")

### set local vol parameters
params <- list(a=0.0012,b=0.1634,sig=0.1029,rho=-0.5555,m1=0.0439)
#params <- list(a=0.04,b=0,sig=0.1029,rho=-0.5555,m1=0.0439)

### local vol MC
strikes <- c(0.8,1.0,1.2) 
print(LocalVolMC(params)(S0=1, capT=1, strikes, N=100000, m=16, evolve=evolveEuler, exactVols=NULL))

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

print(callLocalVolDupirePDE(S0=1, r=0, q=0, sigmaLocalVol, t=1, dK=0.001, dt=0.001, sdw=4, start=1, theta=1/2))