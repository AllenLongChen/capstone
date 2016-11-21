setwd("/home/dongniu/capstone")
source("setup.r")
source("SVI_Fitting/sviFit0.R")
source("SVI_Fitting/sviFitQuarticRoots.r")
source("SVI_Fitting/sviSqrtFit.R")
source("SVI_Fitting/svi.R")
source("LRImpVol2.R")
library("stinepack")

# Now fit SVI-Sqrt
sviSqrtFit <- function(ivolData){
  
  bidVols <- as.numeric(ivolData$Bid);
  askVols <- as.numeric(ivolData$Ask);
  expDates <- unique(ivolData$Texp);
  nSlices <- length(expDates);
  slices <- 1:nSlices;
  sviMatrix <- array(dim=c(nSlices,5));
  #w0 <- computeW0(ivolData);
  nrows <- dim(ivolData)[1];
  midV <- numeric(nrows)*NA;
  kk <- numeric(nrows)*NA;
  ww0 <- numeric(nrows)*NA;
  
  # Compute w0, keeping k and midVar as we go
  for (slice in 1:nSlices){
    t <- expDates[slice];
    texp <- ivolData$Texp;
    bidVol <- bidVols[texp==t];
    askVol <- askVols[texp==t];
    pick <- (bidVol>0)&(askVol>0);
    pick <- !is.na(pick)&(bidVol<askVol);
    midVar <- (bidVol[pick]^2+askVol[pick]^2)/2;
    f <- (ivolData$Fwd[texp==t])[1];
    k <- log(ivolData$Strike[texp==t]/f)[pick]; 
    w0 <- t*stinterp(k,midVar,0)$y;
    # Now put in correct place in columns
    ww0[texp==t][pick]<- w0;
    midV[texp==t][pick]<- midVar;
    kk[texp==t][pick]<- k;    
  }
  
  tcutoff <- min(0.1,max(expDates));
  
  # Define objective function
  obj <- function(sviSqrtParams){
    tmp1 <- as.list(sviSqrtParams);
    names(tmp1) <- c("rho","eta");
    sviSqrtVar <- sviSqrt(tmp1,kk,ww0)/texp;
    tmp <- sum(((midV-sviSqrtVar)^2)[texp >= tcutoff],na.rm=T);
    return(tmp);
  };
  
  sviSqrtGuess <- list(rho = - 0.7, eta = 1.0); 
  
  fit <- optim(sviSqrtGuess,obj,method="L-BFGS-B",lower=c(-.999,-Inf),upper=c(+.999,+Inf));
  res <-as.list(fit$par);
  
  # Now convert the result to an SVI matrix
  sel <- !is.na(ww0);
  w0r <- unique(ww0[sel]);
  rho <- rep(res$rho,nSlices);
  a <- w0r/2*(1-rho^2);
  gg <- res$eta/sqrt(w0r);
  b <- w0r/2*gg;
  m <- -rho/gg;
  sig <- sqrt(1-rho^2)/gg;
  
  tmp <- data.frame(a,b,sig,rho,m);
  
  return(tmp);
  
}



# sviGridToSurface2 <- function(df){
#   df$Fwd <- df$imp_fwd
#   df$Strike <- exp(df$logstrike)*df$Fwd
#   df$Bid <- ifelse(df$logstrike>0,df$call_bid,df$put_bid)
#   df$Ask <- ifelse(df$logstrike>0,df$call_offer,df$put_offer)
#   fitSqrt <- sviSqrtFit(ivolData=df)
#   fitQR <- sviFitQR(ivolData=df, sviGuess=fitSqrt)
#   # print(sviMatrix)
#   impVolSVI <- function(k,t){
#     params <- sviMatrix[sviMatrix$texp==t,c("a","b","sig","rho","m")]
#     # print(params)
#     return(sqrt(impVarSVIRaw(params)(k)/t))
#   }
#   return(impVolSVI)
# }

# sviSurfaces <- function(strikes, texp, trange, df_grid, strike_cutoff = 0.3){
#   
#   svifunc <- sviGridToSurface(subset(df_grid, abs(logstrike) < strike_cutoff))
#   TtlVars <- matrix(, nrow=length(strikes), ncol=1)
#   for(t in texp){
#     
#     temp_iv <- matrix(sapply(strikes, function(x){svifunc(x, t)}), nrow=length(strikes), ncol=1)
#     temp_tv <- temp_iv * temp_iv * t
#     # plot(strikes, temp_iv, xlab="Strikes", ylab="Implied Volatility")
#     # plot(strikes, temp, xlab="Strikes", ylab="Total Variance")
#     TtlVars <- cbind(TtlVars, temp_tv)
#     
#   }
#   TtlVars <- TtlVars[,-1]
#   
#   y <- TtlVars[1, -1]
#   inte <- stinterp(texp[-1], y, trange, method="stineman")
#   time <- inte$x[!is.na(inte$y)]
#   yTV <- inte$y[!is.na(inte$y)]
#   TV_surface <- matrix(yTV)
#   
#   for(i in 2:dim(TtlVars)[1]){
#     
#     y <- TtlVars[i, -1]
#     inte <- stinterp(texp[-1], y, trange, method="stineman")
#     time <- inte$x[!is.na(inte$y)]
#     yTV <- inte$y[!is.na(inte$y)]
#     TV_surface <- cbind(TV_surface, matrix(yTV))
#     
#   }
#   
#   IV_surface <- matrix(, nrow=1, ncol=dim(TV_surface)[2])
#   for(i in 1:dim(TV_surface)[1]){
#     
#     impvols <- matrix(sqrt(TV_surface[i, ] / time[i]), nrow=1, ncol=dim(TV_surface)[2])
#     IV_surface <- rbind(IV_surface, impvols)
#     
#   }
#   IV_surface <- IV_surface[-1, ]
#   
#   persp(time, strikes, TV_surface, col="green", phi=30, theta=150, 
#         r=1/sqrt(3)*20,d=5,expand=.5,ltheta=-135,lphi=20,ticktype="simple",
#         shade=.8,border=NA,ylab="Log-strike k",xlab="Time to Expiration",zlab="Total Variance", main="SVI Fitted Total Variance Surface");
#   
#   persp(time, strikes, IV_surface, col="green", phi=30, theta=30, 
#         r=1/sqrt(3)*20,d=5,expand=.5,ltheta=-135,lphi=20,ticktype="simple",
#         shade=.8,border=NA,ylab="Log-strike k",xlab="Time to Expiration",zlab="Implied Volatility", main="SVI Fitted Implied Volatility Surface");
#   
#   res_surface <- c(TV_surface, IV_surface)
#   return(res_surface)
#   
# }

### subset df_grid by normalized log-strike
strike_subset<-function(df_grid,strike_cutoff=0.9){
  df_grid_subset<-data.frame()
  texp<-unique(df_grid[,"Texp"])
  ### normalized log-strike
  for(t in texp){
    ### atm vol
    df_subset<-subset(df_grid,Texp==t)
    atm_idx<-which(abs(df_subset[,"logstrike"])==min(abs(df_subset[,"logstrike"])))
    atm_vol<-0.5*df_subset[atm_idx,"call_bid"]+0.5*df_subset[atm_idx,"call_offer"]
    ### cutoff by normalized log-strike
    cutoff<-0.9*atm_vol*sqrt(t)
    df_subset<-subset(df_subset,abs(logstrike)<=cutoff)
    df_grid_subset<-rbind(df_grid_subset,df_subset)
  }
  return(df_grid_subset)
}



# generate svi function: only works for existing texp
# uncomment to run
source("/home/dongniu/capstone/SVI_Fitting/sviRoots.R")
source("/home/dongniu/capstone/SVI_Fitting/sviJW.R")
source("/home/dongniu/capstone/SVI_Fitting/sviarbitrage.R")
source("/home/dongniu/capstone/SVI_Fitting/sviVolSurface.R")
source("/home/dongniu/capstone/SVI_Fitting/svi.R")


load("/home/dongniu/capstone/GOOG.RData")
df_grid <- file_to_grid(df)

load("/home/dongniu/capstone/impVols.RData")
df_grid <- impVols
  

df_grid$Fwd <- df_grid$imp_fwd
df_grid$Strike <- exp(df_grid$logstrike)*df_grid$Fwd
df_grid$Bid <- ifelse(df_grid$logstrike>0,df_grid$call_bid,df_grid$put_bid)
df_grid$Ask <- ifelse(df_grid$logstrike>0,df_grid$call_offer,df_grid$put_offer)
df_grid$CallMid <- (df_grid$call_bid + df_grid$call_offer) / 2

df_grid <- strike_subset(df_grid)



fitSqrt <- sviSqrtFit(ivolData=df_grid)
fitQR <- sviFitQR(ivolData=df_grid, sviGuess=fitSqrt, penaltyFactor = 1000)

texp <- df$Texp[!duplicated(df$Texp)]
volTVS <- function(k,t){sqrt(sviW(fitQR,texp,k,t)/t)}


k2 <- seq(-0.3, 0.3, by=0.01) 
t2 <- seq(0.01, 1.4, by=0.02)

z <- t(volTVS(k2, t2))

var <- t(sviW(fitQR, texp, k2, t2))

persp(k2, t2, z, col="green", phi=30, theta=30, 
      r=1/sqrt(3)*20,d=5,expand=.5,ltheta=-135,lphi=20,ticktype="simple",
      shade=.8,border=NA,xlab="Log-strike k",ylab="Expiration t",zlab="Implied vol.");

persp(k2, t2, var, col="green", phi=30, theta=30, 
      r=1/sqrt(3)*20,d=5,expand=.5,ltheta=-135,lphi=20,ticktype="simple",
      shade=.8,border=NA,xlab="Log-strike k",ylab="Expiration t",zlab="Implied vol.");

source("/home/dongniu/capstone/SVI_Fitting/plotIvols.R")
plotIvols(df_grid, fitQR)
plot(t2, var[35,])

