source('LRBinomial.R')
source('BlackScholes.R')

### fit implied vols of selected strikes,
### also get interest rate, borrow cost, implied spot for given Texp
fit_imp_vol0 <- function(df,texp,nStrikes=5,delta=0.2,S0=869.81,sig0=0.2,r0=0,q0=0){
  K_max<-S0*exp((r0-q0+sig0^2/2)*texp-sig0*sqrt(texp)*qnorm(exp((q0)*texp)*delta))
  K_min<-S0*exp((r0-q0+sig0^2/2)*texp+sig0*sqrt(texp)*qnorm(exp((q0)*texp)*delta))
  
  df_subset_c <- subset(df, Texp==texp & cp_flag=="C")
  df_subset_p <- subset(df, Texp==texp & cp_flag=="P")
  
  ### find strike range
  ind_min<- sum(df_subset_c$strike_price<K_min)+1
  ind_max<-sum(df_subset_c$strike_price<K_max)
  
  indx<-seq(ind_min,ind_max,by=floor((ind_max-ind_min)/(nStrikes-1)))
  
  c_px<-(df_subset_c[indx,"best_bid"]+df_subset_c[indx,"best_offer"])/2
  p_px<-(df_subset_p[indx,"best_bid"]+df_subset_p[indx,"best_offer"])/2
  
  c_strike<-df_subset_c[indx,"strike_price"]
  p_strike<-df_subset_p[indx,"strike_price"]
  
  #d1=log(S0*exp((r0-q0)*texp)/c_strike)/sig0/sqrt(texp)+0.5*sig0*sqrt(texp)
  #vega=S0*exp((-q0)*texp)*dnorm(d1)*sqrt(texp)
  
  optim_func<-function(...){
    input<-c(...)
    sig<-input[1:length(c_strike)]
    S<-input[length(c_strike)+1]
    r<-input[length(c_strike)+2]
    q<-input[length(c_strike)+3]
    c_LR<-binPriceLR(c_strike, S, r, q,sig, texp, "C", earlyExercise=TRUE, steps=100)
    p_LR<-binPriceLR(p_strike, S, r, q,sig, texp, "P", earlyExercise=TRUE, steps=100)
    return(sum((c_LR-c_px)^2)+sum((p_LR-p_px)^2))
  }

  # set upper and lower bound for r and q
  r_lower=-0.10
  r_upper=0.30
  q_lower=-0.10
  q_upper=0.30
  sig0=rep(sig0,nStrikes)
  
  res <- optim(c(sig0,S0,r0,q0),optim_func,lower=c(rep(0,nStrikes+1),r_lower,q_lower),
               upper=c(rep(Inf,nStrikes+1),r_upper,q_upper),method='L-BFGS-B',
               control=list(maxit=10000))
  return(c(c_strike,res$par)) # res$par contains the nStrikes imp vols and implied spot, interest rate, borrow cost
}

### fit implied vol smile for given Texp, Call/Put, Bid/Offer
fit_imp_vol<-function(df,texp,cp,bid_offer,S,r,q){
  df_subset <- subset(df, Texp==texp & cp_flag==cp)
  
  opt_prices<-as.numeric(dim(df_subset)[1])
  if(bid_offer=="bid"){
    opt_prices <- df_subset$best_bid
  }
  else{
    opt_prices<-df_subset$best_offer
  }
  
  strikes <- df_subset$strike_price
  
  imp_vols<-as.numeric(length(strikes))

  for(i in 1:length(strikes)){
    strike<-strikes[i]
    opt_price<-opt_prices[i]
    optim_func <- function(sig){
      LRprice<-binPriceLR(strike,S,r,q,sig,texp,cp,earlyExercise=T,steps=100)
      return((LRprice-opt_price)^2)
    }
    sig0<-0.2
    sig1<-optim(sig0,optim_func,lower=0,upper=3,method='Brent',control=list(maxit=1000))
    imp_vols[i]<-sig1$par
  }
  return(data.frame(cbind(strikes,imp_vols)))
}

file_to_grid<-function(df){
  Texps<-unique(df[,"Texp"])
  rates<-as.numeric(length(Texps))
  bcosts<-as.numeric(length(Texps))
  impSpots<-as.numeric(length(Texps))
  nStrikes<-5
  
  impVols<-data.frame()
  
  ### Texp
  ### 3   9  16  23  38  66 129 157 220 521
  
  ### save OTM offer implied vol
  for(i in 1:length(Texps)){
    print(i)
    res<-fit_imp_vol0(df,Texps[i])
    #print(res[1:(2*nStrikes)])
    
    impSpots[i]<-res[length(res)-2]
    rates[i]<-res[length(res)-1]
    bcosts[i]<-res[length(res)]
    cat(impSpots[i],rates[i],bcosts[i])
    imp_vol_c_b<-fit_imp_vol(df,Texps[i],"C","bid",impSpots[i],rates[i],bcosts[i])
    imp_vol_c_o<-fit_imp_vol(df,Texps[i],"C","offer",impSpots[i],rates[i],bcosts[i])
    imp_vol_p_b<-fit_imp_vol(df,Texps[i],"P","bid",impSpots[i],rates[i],bcosts[i])
    imp_vol_p_o<-fit_imp_vol(df,Texps[i],"P","offer",impSpots[i],rates[i],bcosts[i])
    
    df_subset_c <- subset(df, Texp==Texps[i] & cp_flag=="C")
    df_subset_p <- subset(df, Texp==Texps[i] & cp_flag=="P")
    
    fwd<-impSpots[i]*exp((rates[i]-bcosts[i])*Texps[i])
    logstrikes<-log(imp_vol_c_b$strikes/fwd)
    
    imp_vols<-data.frame(matrix(nrow=length(logstrikes),ncol=10))
    colnames(imp_vols)<-c("logstrike","Texp","call_bid","call_offer","put_bid","put_offer",
                          "imp_spot","imp_r","imp_q","imp_fwd")
    imp_vols[,"logstrike"]<-logstrikes
    imp_vols[,"Texp"]<-rep(Texps[i],length(logstrikes))
    imp_vols[,"call_bid"]<-imp_vol_c_b$imp_vols
    imp_vols[,"call_offer"]<-imp_vol_c_o$imp_vols
    imp_vols[,"put_bid"]<-imp_vol_p_b$imp_vols
    imp_vols[,"put_offer"]<-imp_vol_p_o$imp_vols
    imp_vols[,"imp_spot"]<-rep(impSpots[i],length(logstrikes))
    imp_vols[,"imp_r"]<-rep(rates[i],length(logstrikes))
    imp_vols[,"imp_q"]<-rep(bcosts[i],length(logstrikes))
    imp_vols[,"imp_fwd"]<-rep(fwd,length(logstrikes))
    impVols<-rbind(impVols,imp_vols)
    
    #filename=paste(paste(Texps[i]*365,"d",sep=""),".RData",sep="")
    #save(imp_vols,file=filename)
    
    # title=paste(paste("Texp=",Texps[i]*365,sep=""),"d",sep="")
    # ymax=max(max(imp_vol_c_o$imp_vols),max(imp_vol_p_o$imp_vols))
    # 
    # plot(logstrikes,imp_vol_c_b$imp_vols,type="l",xlab="log-strike",ylab="impVol",
    #      col='red',ylim=c(0.0,ymax),main=title)
    # lines(logstrikes,imp_vol_c_o$imp_vols,col="blue")
    # lines(logstrikes,imp_vol_p_b$imp_vols,col="orange")
    # lines(logstrikes,imp_vol_p_o$imp_vols,col="green")
    # points(log(res[1:nStrikes]/fwd),res[(nStrikes+1):(2*nStrikes)])
    # 
    # legend("top",c("call bid","call offer","put bid","put offer","benchmark"),lty=1,col=c("red","blue","orange","green","black"))
  }
  return(impVols)
}
