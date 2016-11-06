source("setup.r")
source("SVI_Fitting/sviFit0.R")
source("SVI_Fitting/sviFitQuarticRoots.r")
source("SVI_Fitting/sviSqrtFit.R")
source("SVI_Fitting/svi.R")
source("LRImpVol2.R")

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
    print(t)
    texp <- ivolData$Texp;
    bidVol <- bidVols[texp==t];
    askVol <- askVols[texp==t];
    pick <- (bidVol>0)&(askVol>0);
    pick <- !is.na(pick);
    midVar <- (bidVol[pick]^2+askVol[pick]^2)/2;
    f <- (ivolData$Fwd[texp==t])[1];
    print(f)
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
  print(sviMatrix)
  impVolSVI <- function(k,t){
    params <- sviMatrix[sviMatrix$texp==t,c("a","b","sig","rho","m")]
    print(params)
    return(sqrt(impVarSVIRaw(params)(k)/t))
  }
  return(impVolSVI)
}

# # generate svi function: only works for existing texp
# # uncomment to run
# load("GOOG.RData")
# df_grid <- file_to_grid(subset(df,Texp==521/365 | Texp==66/365))
# svifunc <- sviGridToSurface(subset(df_grid,abs(logstrike)<0.3))

# # Test
# # #texp <- 521/365
# texp <- 66/365
# strikes <- seq(-1,1,by=0.02)
# plot(strikes,sapply(strikes,function(x){svifunc(x,texp)}),ylim=c(0,0.4))
# lines(df_grid[df_grid$Texp==texp,"logstrike"],df_grid[df_grid$Texp==texp,"call_offer"])
# lines(df_grid[df_grid$Texp==texp,"logstrike"],df_grid[df_grid$Texp==texp,"call_bid"])
# lines(df_grid[df_grid$Texp==texp,"logstrike"],df_grid[df_grid$Texp==texp,"put_offer"])
# lines(df_grid[df_grid$Texp==texp,"logstrike"],df_grid[df_grid$Texp==texp,"put_bid"])
