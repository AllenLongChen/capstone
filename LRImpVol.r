source('setup.r')
source('LRBinomial.R')
source('BlackScholes.R')

load('GOOG.RData')

fit_imp_vol <- function(df,cp,texp,S0,nStrikes){
  df_subset <- subset(df, Texp==texp & cp_flag==cp)
  opt_prices <- as.numeric(0.0*df_subset$best_bid+1.0*df_subset$best_offer) # use mid for now
  strikes <- df_subset$strike_price
  ind <- sum(strikes<S0)
  indx <- (ind-round(nStrikes/2)):(ind+round(nStrikes/2))
#   indx <- seq(1,length(strikes),by=round(length(strikes)/nStrikes))
  opt_prices <- opt_prices[indx]
  strikes <- strikes[indx]
  type <- cp_flag
  optim_func <- function(...){
    input <- c(...)
    sig <- input[1:(length(input)-1)]
    r <- input[length(input)]
    LRprice <- function(sig,r){
      return(as.numeric(binPriceLR(strikes,S0,r,sig,texp,type,earlyExercise=T,steps=100)))
    }
    return(sum((LRprice(sig,r)-opt_prices)^2))  
  }
  # initial values
  sig0 <- rep(0.1,length(strikes))
  r0 <- 0.00
  r_upper <- 0.0001
  r_lower <- 0.00
  res <- optim(c(sig0,r0),optim_func,lower=c(rep(-Inf,length(sig0)),r_lower),
               upper=c(rep(Inf,length(sig0)),r_upper),method='L-BFGS-B')
  res$strikes <- strikes
  res$opt_prices <- opt_prices
  return(res) # res$par contains the nStrikes imp vols and the interest rate
}

cp_flag <- 'C'
texp <- 38/365
S0 <- 869.81
nStrikes <- 10
res <- fit_imp_vol(df,cp_flag,texp,S0,nStrikes)
plot(res$strikes,res$par[1:length(res$strikes)],ylim=c(-0.1,0.3))
points(res$strikes,BSImpliedVolCall(S0,res$strikes,texp,0,res$opt_prices),col='red')
# note the function for imp

