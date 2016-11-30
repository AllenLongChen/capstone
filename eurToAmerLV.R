# setwd("/home/dongniu/capstone")
source("BlackScholes.R")
source("SVILocalVolToEurImpVol.R")

#######################################################
# params       is a list that consists of r, S0, iscall
# localVolFunc is a function that takes time to expiry 
#              and strike as arguments and return
#              a single local volatility value
# impVolFunc   is a function that takes time to expiry 
#              and strike as arguments and return
#              a single local volatility value
#######################################################

# formula 3.10
# returns a function that takes valuation time t, 
# expiry capT, and strike K
getEarlyExerciseBoundaryFunc <- function(params,localVolFunc){
  r <- params$r
  ### early exercise boundary for put option
  
  func <- function(t,capT,K){
    tau <- capT - t
    sigl <- localVolFunc(tau,K)
    gamma <- sigl/8/pi/r^2/K^2
    
    res <- K-sigl*sqrt(tau*log(gamma/tau*(1-2/log(tau/gamma))))
    
    res <- max(res,0.0001)
    return(res)
  }
  return(func)
}

# European option price function given implied vol
# returns a function that takes valuation time t
# expiry capT, and strike K
# P(T,K) can be obtained by calling
# getEuropeanPriceFunc(params,impVolFunc)(0,T,K)
getEuropeanPriceFunc <- function(params,impVolFunc){
  r <- params$r
  S0 <- params$S0
  iscall <- params$iscall
  if (iscall) {return( function(t,capT,K){BSFormula(S0,K,capT-t,r,impVolFunc(capT-t,K))} )}
  else {return( function(t,capT,K){BSFormulaPut(S0,K,capT-t,r,impVolFunc(capT-t,K))} )}  
}

# formula 3.11
# returns a function that takes time to expiry and strike
getAmericanPriceFunc <- function(params,localVolFunc,impVolFunc){
  getEarlyExercisePremium <- function(Texp,K){
    S0 <- params$S0
    r  <- params$r
    dPdKFunc <- function(params){
      calc <- function(t,capT,K){
        sigma <- impVolFunc(capT-t,K)
        x <- log(S0/K)+r*(capT-t)
        sig <- sigma*sqrt(capT-t)
        d1 <- x/sig+sig/2
        d2 <- d1 - sig
        #return(pnorm(-d2))
        return(pnorm(-d2)*exp(-r*(capT-t)))
      }
      return(calc)
    }
    dt <- Texp/100
    ############### two 0.01s that need to be taken care of
    ts <- seq(0.0005,Texp,by=dt)
    integrands <- sapply(ts,function(t){dPdKFunc(params)(0,t,getEarlyExerciseBoundaryFunc(params,localVolFunc)(t,Texp,K))})
    ##############
    return(r*K*sum(1/2*(integrands[-1]+integrands[-length(integrands)])*dt))
  }
  
  americanFunc <- function(Texp,K){
    if(!params$iscall){
      euroOptPrice <- getEuropeanPriceFunc(params,impVolFunc)(0,Texp,K)
      earlyExercisePrem <- getEarlyExercisePremium(Texp,K)
      return(euroOptPrice+earlyExercisePrem)
    }
    else{
      return(getEuropeanPriceFunc(params,impVolFunc)(0,Texp,K))
    }
    
  }
  return(americanFunc)
}

# Test region

#params <- list(S0 = 1, r = 0.02, iscall = FALSE)
#earlybound <- getEarlyExerciseBoundaryFunc(params, SVI_LocalVol_func(svi_param, params$S0))
#eurprice <- getEuropeanPriceFunc(params, EurImpVol_func(svi_param, params$S0))
#ameprice <- getAmericanPriceFunc(params, SVI_LocalVol_func(svi_param, params$S0), EurImpVol_func(svi_param, params$S0))


