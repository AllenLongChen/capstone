load("GOOG.RData")
put<-subset(df,cp_flag=="P")
put<-subset(df,Texp==521/365)
put[,"mid"]<-(put[,"best_bid"]+put[,"best_offer"])/2
put<-put[,c("strike_price","Texp","mid")]

source("eurToAmerLV.R")

optimize<-function(df){
  optim_func<-function(...){
    input <- c(...)
    svi_param<-as.list(input)
    names(svi_param)<-c("a","b","sigma","rho","m")
    params <- list(S0 = 869.81, r = 0.01, iscall = FALSE)
    ameprice <- getAmericanPriceFunc(params, SVI_LocalVol_func(svi_param,params$S0), EurImpVol_func(svi_param, params$S0))
    ameprice1<-apply(put[,c("Texp","strike_price")],1,function(y){ameprice(y["Texp"],y["strike_price"])})
    return(sum((ameprice1-put[,"mid"])^2))
  }
  svi_param0<-c(0.0012, 0.1634,0.1029,0.5555,0.0439)
  res <- optim(svi_param0,optim_func,control=list(maxit=10000))
  return(res$par)
}

svi_param_fitted<-optimize(put)
print(svi_param_fitted)