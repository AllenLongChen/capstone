BSFormula <- function(S0, K, capT, r, sigma)
{
x <- log(S0/K)+r*capT;
sig <- sigma*sqrt(capT);
d1 <- x/sig+sig/2;
d2 <- d1 - sig;
pv <- exp(-r*capT);
return( S0*pnorm(d1) - pv*K*pnorm(d2));
}

BSFormulaPut <- function(S0, K, capT, r, sigma)
{
x <- log(S0/K)+r*capT;
sig <- sigma*sqrt(capT);
d1 <- x/sig+sig/2;
d2 <- d1 - sig;
pv <- exp(-r*capT);
return( S0*pnorm(d1) - pv*K*pnorm(d2)+pv*K-S0);
}

# This function now works with vectors of strikes and option values
BSImpliedVolCall <- function(S0, K, capT, r, C)
{
nK <- length(K);
sigmaL <- rep(1e-10,nK);
CL <- BSFormula(S0, K, capT, r, sigmaL);
sigmaH <- rep(10,nK);
CH <- BSFormula(S0, K, capT, r, sigmaH);

  while (mean(sigmaH - sigmaL) > 1e-10)
  {
    sigma <- (sigmaL + sigmaH)/2;
    CM <- BSFormula(S0, K, capT, r, sigma);
    CL <- CL + (CM < C)*(CM-CL);
    sigmaL <- sigmaL + (CM < C)*(sigma-sigmaL);
    CH <- CH + (CM >= C)*(CM-CH);
    sigmaH <- sigmaH + (CM >= C)*(sigma-sigmaH);
  }
  return(sigma);
}

# This function also works with vectors of strikes and option values  
BSImpliedVolPut <- function(S0, K, capT, r, P)
{
#pv <- exp(-r*capT);
sigmaL <- 1e-10;
#intrinsic <- (K-pv*S0);
nK <- length(K);
sigmaL <- rep(1e-10,nK);
#PL <- BSFormula(S0, K, capT, r, sigmaL)+intrinsic;
PL <- BSFormula(S0, K, capT, r, sigmaL);
sigmaH <- rep(10,nK);
#PH <- BSFormula(S0, K, capT, r, sigmaH)+intrinsic;
PH <- BSFormula(S0, K, capT, r, sigmaH);
  while (mean(sigmaH - sigmaL) > 1e-10)
  {
    sigma <- (sigmaL + sigmaH)/2;
    #PM <- BSFormula(S0, K, capT, r, sigma)+intrinsic;
    PM <- BSFormulaPut(S0, K, capT, r, sigma);
    PL <- PL + (PM < P)*(PM-PL);
    sigmaL <- sigmaL + (PM < P)*(sigma-sigmaL);
    PH <- PH + (PM >= P)*(PM-PH);
    sigmaH <- sigmaH + (PM >= P)*(sigma-sigmaH);
  }
  return(sigma);
}

# Function to compute option prices and implied vols given list of final values of underlying
bsOut <- function(xf,capT,AK)
{
    nK <- length(AK);
    N <- length(xf);
    xfbar <- mean(xf);
    CAV <- numeric(nK); BSV <- numeric(nK);
    BSVL <- numeric(nK); BSVH <- numeric(nK);
    for (j in 1:nK){
                    payoff <- (xf-AK[j]) * (xf>AK[j]);
                    CAV[j] <- sum(payoff)/N;
                    err <- sqrt(var(payoff)/N);
                    BSV[j] <- BSImpliedVolCall(xfbar, AK[j], capT,0, CAV[j]);
                    BSVL[j] <- BSImpliedVolCall(xfbar, AK[j], capT,0, CAV[j]-err);
                    BSVH[j] <- BSImpliedVolCall(xfbar, AK[j], capT,0, CAV[j]+err);
                    }
    return(data.frame(AK,CAV,BSV,BSVL,BSVH));
    
}

# Function to return implied vols for a range of strikes
analyticOut <- function(callFormula,AK,capT)
{
    nK <- length(AK);
    #callFormula is a function that computes the call price
    callPrice <- numeric(nK); BSV <- numeric(nK);
    for (j in 1:nK){
                    callPrice[j] <- callFormula(AK[j]);
                    BSV[j] <- BSImpliedVolCall(1, AK[j], capT,0, callPrice[j]);
                    }
    return(data.frame(AK,callPrice,BSV));
}
