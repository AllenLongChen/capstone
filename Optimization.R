setwd("/Users/AllenLongChen/git/capstone/")

load("GOOG.RData")
put<-subset(df,cp_flag=="P")
put<-subset(put,Texp==521/365)
put[,"mid"]<-(put[,"best_bid"]+put[,"best_offer"])/2
put<-put[,c("strike_price","Texp","mid")]
put[,"mid"]<-put[,"mid"]/869.81 
put[,"strike_price"]<-put[,"strike_price"]/869.81

source("eurToAmerLV.R")

optimize<-function(df){
  optim_func<-function(...){
    input <- c(...)
    svi_param<-as.list(input)
    print(svi_param)
    names(svi_param)<-c("a","b","sigma","rho","m")
    params <- list(S0 = 1, r = 0.01, iscall = FALSE)
    ameprice <- getAmericanPriceFunc(params, SVI_LocalVol_func(svi_param,params$S0), EurImpVol_func(svi_param, params$S0))
    ameprice1<-apply(put[,c("Texp","strike_price")],1,function(y){ameprice(y["Texp"],y["strike_price"])})
    return(sum((ameprice1-put[,"mid"])^2))
  }
  svi_param0<-c(0.0012, 0.1634,0.1029,-0.5555,0.0439)
  res <- optim(svi_param0,optim_func,control=list(maxit=10000),method="L-BFGS-B",lower=c(-Inf,0.0001,-Inf,-Inf,-Inf))
  return(res$par)
}

svi_param_fitted<-optimize(put)
param_names <- c("a","b","sigma","rho","m")
svi_param_fitted_list <- list()
for (i in 1:5){
  svi_param_fitted_list[param_names[i]] <- svi_param_fitted[i]
}
print(svi_param_fitted_list)
