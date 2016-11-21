# SVI parameterization of the local volatility surface

svi_param <- list(a = 0.0012, b = 0.1634, sigma = 0.1029, rho = -0.5555, m = 0.0439)

SVI_LocalVol_slice <- function(param, k, texp){
  a <- param$a
  b <- param$b
  sigma <- param$sigma
  rho <- param$rho
  m <- param$m

  LocVol2 <- a + b * (rho * (k / sqrt(texp) - m) + sqrt((k / sqrt(texp) - m)^2 + sigma^2 * texp))

  return(LocVol2)
}

SVI_LocalVol_surface <- function(param, k, time){
  
  surface <- sapply(time, function(x){sqrt(SVI_LocalVol_slice(param, k, x))})
  return(surface)
}

SVI_LocalVol_func <- function(param, X0){
  
  func <- function(time, K){
    k <- log(K / X0)
    surface <- sapply(time, function(x){sqrt(SVI_LocalVol_slice(param, k, x))})
    return(surface)
  }
  return(func)
}

# SVI_LocalVol_func <- function(param){
#   
#   a <- param$a
#   b <- param$b
#   sigma <- param$sigma
#   rho <- param$rho
#   m <- param$m
#   func <- function(texp,k){
#     LocVol2 <- a + b * (rho * (k / sqrt(texp) - m) + sqrt((k / sqrt(texp) - m)^2 + sigma^2 * texp))
#     return(LocVol2)
#   }
#   
#   return(func)
#   
# }


# European implied volatility from Dupire's estimation (3.6)
# param encloses all the SVI local vol parameters
# X0 is the spot price at time 0
# time range and k = log(K/X0) are needed for specific vols

EurImpVol_func <- function(param, X0){
  EurImpVol_slice1 <- function(param, k, texp){
    ImpVol <- c()
    for(i in 1:length(k)){
      # ATM implied vol equals to ATM local vol, no calculation required
      if(k[i] == 0){
        ImpVol <- c(ImpVol, sqrt(SVI_LocalVol_slice(param, 0, texp)))
        next
      }
      
      K <- exp(k[i]) * X0
      inverse <- function(x){
        1 / sqrt(SVI_LocalVol_slice(param, log(x / X0), texp))
      }
      denominator <- integrate(inverse, X0, K, stop.on.error = FALSE)$value
      ImpVol <- c(ImpVol, k[i] / denominator)
    }
    return(ImpVol)
  }
  func <- function(time, K){
    k <- log(K / X0)
    surface <- sapply(time, function(x){EurImpVol_slice1(param, k, x)})
    return(surface)
  }
  return(func)
}


# EurImpVol_surface <- function(param, k, X0, time){
#   
#   EurImpVol_slice <- function(param, k, X0, texp){
#     
#     ImpVol <- c()
#     
#     for(i in 1:length(k)){
#       
#       if(k[i] == 0){
#         
#         ImpVol <- c(ImpVol, sqrt(SVI_LocalVol_slice(param, 0, texp)))
#         next
#         
#       }
#       
#       K <- exp(k[i]) * X0
#       inverse <- function(x){
#         
#         1 / sqrt(SVI_LocalVol_slice(param, log(x / X0), texp))
#         
#       }
#       
#       denominator <- integrate(inverse, X0, K)$value
#       ImpVol <- c(ImpVol, k[i] / denominator)
#       
#     }
#     return(ImpVol)
#     
#   }
#   
#   surface <- sapply(time, function(x){EurImpVol_slice(param, k, X0, x)})
#   return(surface)
#   
# }
# 


# Test calculation region
# 
# k <- seq(-1.0, 1.0, by = 0.02)
# time <- seq(0.1, 1.5, by = 0.1)
# 
# localvol <- SVI_LocalVol_surface(param, k, time)
# impvol <- EurImpVol_surface(param, k, 1, time)
# 
# persp(k, time, localvol, col="green", phi=30, theta=30, 
#       r=1/sqrt(3)*20,d=5,expand=.5,ltheta=-135,lphi=20,ticktype="simple",
#       shade=.8,border=NA,,ylab="Time to Expiration",xlab="Log-Strike k",zlab="Local Volatility", main="SVI Local Vol Surface")
# 
# persp(k, time, impvol, col="green", phi=30, theta=30, 
#       r=1/sqrt(3)*20,d=5,expand=.5,ltheta=-135,lphi=20,ticktype="simple",
#       shade=.8,border=NA,,ylab="Time to Expiration",xlab="Log-Strike k",zlab="Implied Volatility", main="Implied Vol Surface")
