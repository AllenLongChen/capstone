source("BlackScholes.R")

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
getEarlyExerciseBoundaryFunc(params,localVolFunc){
  r <- params$r
  func <- function(t,capT,K){
    tau <- capT - t
    sigl <- localVolFunc(tau,K)
    gamma <- sigl/8/pi/r^2/K^2
    return( K-sigl*sqrt(tau*log(gamma/tau*(1-2/log(tau/gamma)))) )
  }
  return(func)
}

# European option price function given implied vol
# returns a function that takes valuation time t
# expiry capT, and strike K
# P(T,K) can be obtained by calling
# getEuropeanPriceFunc(params,impVolFunc)(0,T,K)
getEuropeanPriceFunc(params,impVolFunc){
  r <- params$r
  S0 <- params$S0
  iscall <- params$iscall
  if (iscall) {return( function(t,capT,K){BSFormula(S0,K,capT-t,r,impVolFunc(capT-t,K))} )}
  else {return( function(t,capT,K){BSFormulaPut(S0,K,capT-t,r,impVolFunc(capT-t,K))} )}  
}

# formula 3.11
# returns a function that takes time to expiry and strike
getAmericanPriceFunc(params,localVolFunc,impVolFunc){
  getEarlyExercisePremium <- function(Texp,K){
    S0 <- params$S0
    r  <- params$r
    dPdKFunc <- function(params){
      calc <- function(t,capT,K){
        sig <- impVolFunc(capT-t,K)
        x <- log(S0/K)+r*(capT-t)
        sig <- sigma*sqrt(capT-t)
        d1 <- x/sig+sig/2
        d2 <- d1 - sig
        return(pnorm(-d2))
      }
      return(calc)
    }
    dt <- Texp/100
    ts <- seq(0,Texp,by=dt)
    integrands <- sapply(ts,function(t){dPdKFunc(params)(0,t,getEarlyExerciseBoundaryFunc(params,localVolFunc)(t,Texp,K))})
    return(r*K*sum(1/2*(integrands[-1]+integrands[-length(integrands)])*dt))
  }
  
  americanFunc <- function(Texp,K){
    euroOptPrice <- getEuropeanPriceFunc(params,impVolFunc)(0,Texp,K)
    earlyExercisePrem <- getEarlyExercisePremium(Texp,K)
    return(euroOptPrice+earlyExercisePrem)
  }
  return(americanFunc)
}
