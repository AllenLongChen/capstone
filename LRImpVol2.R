source('LRBinomial.R')
source('BlackScholes.R')

load('GOOG.RData')

### fit implied vols of selected strikes,
### also get interest rate, borrow cost, implied spot for given Texp
fit_imp_vol0 <- function(df,texp,nStrikes=5,delta=0.2,S0=869.81,sig0=0.2,r0=0,q0=0){
  K_max<-S0*exp((r0-q0+sig0^2/2)*texp-sig0*sqrt(texp)*qnorm(exp(q0*texp)*delta))
  K_min<-S0*exp((r0-q0+sig0^2/2)*texp+sig0*sqrt(texp)*qnorm(exp(q0*texp)*delta))
  
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
  r_lower=0.00
  r_upper=0.20
  q_lower=0.00
  q_upper=0.20
  sig0=rep(sig0,nStrikes)
  
  res <- optim(c(sig0,S0,r0,q0),optim_func,lower=c(rep(0,nStrikes+1),r_lower,q_lower),
               upper=c(rep(Inf,nStrikes+1),r_upper,q_upper),method='L-BFGS-B',
               control=list(maxit=10000))
  return(res) # res$par contains the nStrikes imp vols and implied spot, interest rate, borrow cost
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
    sig1<-optim(sig0,optim_func,lower=0.1,upper=0.9,method='L-BFGS-B',control=list(maxit=1000))
    imp_vols[i]<-sig1$par
  }
  #print(imp_vols)
  return(data.frame(cbind(strikes,imp_vols)))
}

Texps<-unique(df[,"Texp"])
rates<-as.numeric(length(Texps))
bcosts<-as.numeric(length(Texps))
impSpots<-as.numeric(length(Texps))

### Texp
### 3   9  16  23  38  66 129 157 220 521

for(i in 1:length(Texps)){
  if(i==10){
    res<-fit_imp_vol0(df,Texps[i])$par
    print(res)
    impSpots[i]<-res[length(res)-2]
    rates[i]<-res[length(res)-1]
    bcosts[i]<-res[length(res)]
    imp_vol_c_b<-fit_imp_vol(df,Texps[i],"C","bid",impSpots[i],rates[i],bcosts[i])
    imp_vol_c_o<-fit_imp_vol(df,Texps[i],"C","offer",impSpots[i],rates[i],bcosts[i])
    imp_vol_p_b<-fit_imp_vol(df,Texps[i],"P","bid",impSpots[i],rates[i],bcosts[i])
    imp_vol_p_o<-fit_imp_vol(df,Texps[i],"P","offer",impSpots[i],rates[i],bcosts[i])
    plot(imp_vol_c_b$strikes,imp_vol_c_b$imp_vols,type="l",xlab="strike",ylab="impVol",
         col='red',ylim=c(0.14,0.4))
    lines(imp_vol_c_o$strikes,imp_vol_c_o$imp_vols,col="blue")
    lines(imp_vol_p_b$strikes,imp_vol_p_b$imp_vols,col="orange")
    lines(imp_vol_p_o$strikes,imp_vol_p_o$imp_vols,col="green")
    legend("topright",c("call bid","call offer","put bid","put offer"),lty=1,col=c("red","blue","orange","green"))
  }
  # res<-fit_imp_vol0(df,Texps[i])$par
  # impSpots[i]<-res[length(nStrikes)+1]
  # rates[i]<-res[length(nStrikes)+2]
  # bcosts[i]<-res[length(nStrieks)+3]
}