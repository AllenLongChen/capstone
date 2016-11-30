load("GOOG.RData")

df[,"mid"]<-(df[,"best_bid"]+df[,"best_offer"])/2

df<-df[,c("strike_price","Texp","mid","cp_flag")]
df[,"mid"]<-df[,"mid"]/869.81 
df[,"strike_price"]<-df[,"strike_price"]/869.81

df_c<-subset(df,cp_flag=="C")

df_p<-subset(df,cp_flag=="P")

df_c<-subset(df_c,strike_price>1)
df_p<-subset(df_p,strike_price<1)

# subindex_p<-seq(5,nrow(df_p),10)
# subindex_c<-seq(5,nrow(df_c),10)
# df_c<-df_c[subindex_c,]
# df_p<-df_p[subindex_p,]

source("eurToAmerLV.R")

Texps<-c(3,9,16,23,38,66,129,157,220,521)/365

optimize<-function(df_c,df_p){
  optim_func<-function(...){
    input <- c(...)
    print(input)
    svi_param<-as.list(input[1:5])
    names(svi_param)<-c("a","b","sigma","rho","m")
    rs<-input[6:15]
    
    ### add penalty to negative vols
    penalty<-svi_param$a+svi_param$b*svi_param$sigma*sqrt(1-svi_param$rho^2)*unique(df[,"Texp"])
    print(min(penalty))
    if(min(penalty)<0){
      return(1000000*(-min(penalty)))
    }
    
    ### calculate cost function
    res<-0
    
    for(i in 1:10){
      ### put
      params_put <- list(S0 = 1, r = rs[i], iscall = FALSE)
      ### call
      params_call<-list(S0=1,r=rs[i],iscall=TRUE)
      
      df_p1<-subset(df_p,Texp==Texps[i])
      df_c1<-subset(df_c,Texp==Texps[i])
      
      ameprice_put <- getAmericanPriceFunc(params_put, SVI_LocalVol_func(svi_param,params_put$S0), EurImpVol_func(svi_param, params_put$S0))
      ameprice1_put<-sapply(df_p1["strike_price"],function(y){ameprice_put(Texps[i],y)})
      
      ameprice_call <- getAmericanPriceFunc(params_call, SVI_LocalVol_func(svi_param,params_call$S0), EurImpVol_func(svi_param, params_call$S0))
      ameprice1_call<-sapply(df_c1["strike_price"],function(y){ameprice_call(Texps[i],y)})
      
      res<-res+sum((ameprice1_put-df_p1[,"mid"])^2)+sum((ameprice1_call-df_c1[,"mid"])^2)
    }
    return(res)
    
  }
  svi_param0<-c(0.0012, 0.1634,0.1029,-0.5555,0.0439,rep(0.01,10))
  res <- optim(svi_param0,optim_func,method="L-BFGS-B",lower=c(-Inf,1e-9,1e-9,-1,-Inf,rep(-0.1,10)),upper=c(Inf,Inf,Inf,1,Inf,rep(0.1,10)))
  
  return(res$par)
}

svi_param_fitted<-optimize(df_c,df_p)
param_names <- c("a","b","sigma","rho","m")
svi_param_fitted_list <- list()
for (i in 1:5){
  svi_param_fitted_list[param_names[i]] <- svi_param_fitted[i]
}
print(svi_param_fitted_list)
