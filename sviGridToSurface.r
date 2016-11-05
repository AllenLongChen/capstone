setwd("~/git/capstone/")
source("setup.r")
source("SVI_Fitting/sviFit0.R")
source("SVI_Fitting/sviFitQuarticRoots.r")
source("SVI_Fitting/sviSqrtFit.R")
source("SVI_Fitting/svi.R")
source("LRImpVol2.R")

impVolSVIRaw <- function(params){
  return(function(k){svi(params,k)})
}

sviGridToSurface <- function(df){
#   texps <- unique(df$Texp)
#   logstrike
#   imp_fwd
#   call_bid
#   call_offer
#   put_bid
#   put_offer
  df$Bid <- ifelse(df$logstrike>0,df$call_bid,df$put_bid)
  df$Ask <- ifelse(df$logstrike>0,df$call_offer,df$put_offer)
  sviMatrix <- sviFit(ivolData=df)
  
#   applySVI <- function(df){
#     'svi function'
#     'list=assign params  '
#     'params$texp <- '
#     'params$a <- '
#     'params$b <- '
#     'params$rho <- '
#     'params$m <- '
#     'params$sigma <- '
#     return(params)
#   }
#   res <- data.frame(texp=character(0),a=numeric(0),b=numeric(0),rho=numeric(0),m=numeric(0),sigma=numeric(0))
#   for (texp in texps){
#     df_subset <- subset(df, 'time to expiry'==texp)
#     res <- rbind(res,applySVI(df_subset))
#   }
# #   return(res)
  
#   impVolSVI <- function(k,t){
#     params <- 'from res find texp=t'
#     'sviFormula(params)'
#   }
#   return(impVolSVI)
  return(sviMatrix)
}

load("GOOG.RData")
df_grid <- file_to_grid(subset(df,Texp==521/365 & abs(strike_price-870)<200))
df_params <- sviGridToSurface(df_grid)
strikes <- seq(-1,1,by=0.02)
plot(strikes,impVolSVIRaw(as.list(df_params))(strikes),ylim=c(0,0.4))
lines(df_grid$logstrike,df_grid$call_offer)
lines(df_grid$logstrike,df_grid$call_bid)
lines(df_grid$logstrike,df_grid$put_bid)
lines(df_grid$logstrike,df_grid$put_offer)
