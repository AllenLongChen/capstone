setwd("/home/dongniu/capstone")
source("setup.r")
source("SVI_Fitting/sviFit0.R")
source("SVI_Fitting/sviFitQuarticRoots.r")
source("SVI_Fitting/sviSqrtFit.R")
source("SVI_Fitting/svi.R")
source("LRImpVol2.R")
library("stinepack")

impVarSVIRaw <- function(params){
  return(function(k){svi(params,k)})
}

sviFit <- function(ivolData){
  
  bidVols <- as.numeric(ivolData$Bid);
  askVols <- as.numeric(ivolData$Ask);
  expDates <- unique(ivolData$Texp);
  nSlices <- length(expDates);
  slices <- 1:nSlices;
  sviMatrix <- array(dim=c(nSlices,6));
  
  ######################################
  for (slice in slices){
    t <- expDates[slice];
   # print(t)
    texp <- ivolData$Texp;
    bidVol <- bidVols[texp==t];
    askVol <- askVols[texp==t];
    pick <- (bidVol>0)&(askVol>0);
    pick <- !is.na(pick);
    midVar <- (bidVol[pick]^2+askVol[pick]^2)/2;
    f <- (ivolData$Fwd[texp==t])[1];
   # print(f)
    k <- log(ivolData$Strike[texp==t]/f)[pick]; # Plot vs log-strike
    plot(k,bidVol);
    lines(k,askVol);
    
    sviGuess <- list(a = mean(midVar,na.rm=T), b = 0.1, sig = 0.1, rho = - 0.7, m = 0); 
    
    # Define objective function
    obj <- function(sviparams){
      tmp1 <- as.list(sviparams);
      names(tmp1) <- c("a","b","sig","rho","m")
      sviVar <- svi(tmp1,k);
      minVar <- tmp1$a+tmp1$b*tmp1$sig*sqrt(abs(1-tmp1$rho^2));
      negVarPenalty <- min(100,exp(-1/minVar));
      tmp <- sum((midVar-sviVar)^2,na.rm=T)+negVarPenalty;
      return(tmp*10000);
    };
    
    fit <- optim(sviGuess,obj);
    if(abs(fit$par["rho"])>.999){
      fit <- optim(fit$par,obj,method="L-BFGS-B",lower=c(-10,0,0.00000001,-.999,-10),upper=c(+10,100,100,+.999,+10));}
    sviMatrix[slice,] <- c(t,fit$par*c(t,t,1,1,1));
    
  }# end of for loop
  ######################################
  sviMatrix <- as.data.frame(sviMatrix);
  colnames(sviMatrix) <- c("texp","a","b","sig","rho","m");
  return(sviMatrix);
}

sviGridToSurface <- function(df){
  df$Fwd <- df$imp_fwd
  df$Strike <- exp(df$logstrike)*df$Fwd
  df$Bid <- ifelse(df$logstrike>0,df$call_bid,df$put_bid)
  df$Ask <- ifelse(df$logstrike>0,df$call_offer,df$put_offer)
  sviMatrix <- sviFit(ivolData=df)  
  # print(sviMatrix)
  impVolSVI <- function(k,t){
    params <- sviMatrix[sviMatrix$texp==t,c("a","b","sig","rho","m")]
   # print(params)
    return(sqrt(impVarSVIRaw(params)(k)/t))
  }
  return(impVolSVI)
}

sviSurfaces <- function(strikes, texp, trange, df_grid, strike_cutoff = 0.3){
  
  svifunc <- sviGridToSurface(subset(df_grid, abs(logstrike) < strike_cutoff))
  TtlVars <- matrix(, nrow=length(strikes), ncol=1)
  for(t in texp){
    
    temp_iv <- matrix(sapply(strikes, function(x){svifunc(x, t)}), nrow=length(strikes), ncol=1)
    temp_tv <- temp_iv * temp_iv * t
    # plot(strikes, temp_iv, xlab="Strikes", ylab="Implied Volatility")
    # plot(strikes, temp, xlab="Strikes", ylab="Total Variance")
    TtlVars <- cbind(TtlVars, temp_tv)
    
  }
  TtlVars <- TtlVars[,-1]
  
  y <- TtlVars[1, -1]
  inte <- stinterp(texp[-1], y, trange, method="stineman")
  time <- inte$x[!is.na(inte$y)]
  yTV <- inte$y[!is.na(inte$y)]
  TV_surface <- matrix(yTV)
  
  for(i in 2:dim(TtlVars)[1]){
    
    y <- TtlVars[i, -1]
    inte <- stinterp(texp[-1], y, trange, method="stineman")
    time <- inte$x[!is.na(inte$y)]
    yTV <- inte$y[!is.na(inte$y)]
    TV_surface <- cbind(TV_surface, matrix(yTV))
    
  }
  
  IV_surface <- matrix(, nrow=1, ncol=dim(TV_surface)[2])
  for(i in 1:dim(TV_surface)[1]){
    
    impvols <- matrix(sqrt(TV_surface[i, ] / time[i]), nrow=1, ncol=dim(TV_surface)[2])
    IV_surface <- rbind(IV_surface, impvols)
    
  }
  IV_surface <- IV_surface[-1, ]
  
  persp(time, strikes, TV_surface, col="green", phi=30, theta=150, 
        r=1/sqrt(3)*20,d=5,expand=.5,ltheta=-135,lphi=20,ticktype="simple",
        shade=.8,border=NA,ylab="Log-strike k",xlab="Time to Expiration",zlab="Total Variance", main="SVI Fitted Total Variance Surface");
  
  persp(time, strikes, IV_surface, col="green", phi=30, theta=30, 
        r=1/sqrt(3)*20,d=5,expand=.5,ltheta=-135,lphi=20,ticktype="simple",
        shade=.8,border=NA,ylab="Log-strike k",xlab="Time to Expiration",zlab="Implied Volatility", main="SVI Fitted Implied Volatility Surface");
  
  res_surface <- c(TV_surface, IV_surface)
  return(res_surface)
  
}

 # generate svi function: only works for existing texp
 # uncomment to run
 load("/home/dongniu/capstone/GOOG.RData")
 df_grid <- file_to_grid(df)
 # svifunc <- sviGridToSurface(subset(df_grid,abs(logstrike)<0.3))
 
 texp <- df$Texp[!duplicated(df$Texp)]
 strikes <- seq(-1,1,by=0.02) 
 trange <- seq(0, 10, by=0.01)
 
 TestSurface <- sviSurfaces(strikes, texp, trange, df_grid)
 
 # Test
 # #texp <- 521/365
 # texp1 <- 66/365
 # texp2 <- 521/365
 # p1 <- sapply(strikes, function(x){svifunc(x,texp1)})
 # p2 <- sapply(strikes, function(x){svifunc(x,texp2)})
 # plot(strikes, p1, ylim=c(0, 0.6))
 # plot(strikes, p2, ylim=c(0, 0.6))
 # x <- c(texp1, texp2)
 # y <- c(p1[52], p2[52])

 # lines(df_grid[df_grid$Texp==texp,"logstrike"],df_grid[df_grid$Texp==texp,"call_offer"])
 # lines(df_grid[df_grid$Texp==texp,"logstrike"],df_grid[df_grid$Texp==texp,"call_bid"])
 # lines(df_grid[df_grid$Texp==texp,"logstrike"],df_grid[df_grid$Texp==texp,"put_offer"])
 # lines(df_grid[df_grid$Texp==texp,"logstrike"],df_grid[df_grid$Texp==texp,"put_bid"])
