fh <- function(z, n)
{
  return (0.5 + sign(z)*sqrt(0.25-0.25*exp(-((z/(n+1/3+0.1/(n+1)))^2)*(n+1/6))))
}

### modified to include borrow cost
binPriceLR <- function (X, S0, r, q,sig, Texp, oType, earlyExercise=TRUE, steps=10000)
{
  # Function to calculate the price of a vanilla European or American
  # Put or Call option using a Leisen-Reimer binomial tree.
  #
  # Inputs: X - strike
  #       : S0 - stock price
  #       : r - risk free interest rate
  #       : q - borrow cost
  #       : sig - volatility
  #       : Texp -  Time to expiry
  #       : steps - number of time steps to calculate
  #       : oType - must be 'PUT' or 'CALL'.
  #       : earlyExercise - true for American, false for European.
  #
  # Output: oPrice - the option price
  #
  
  # Calculate the d1 and d2 parameters of the Black-Scholes model
  #T <- dt*steps # time to maturity
  
  nK <- length(X)
  capT <- Texp
  dt <- Texp / steps
  
  # Calculate the Leisen-Reimer model parameters
  # The parameter n must be odd
  n <- 0
  if((steps %% 2) > 0)
  {
    n <- steps
  }
  else
  {
    n <- steps+1
  }
  
  V <- numeric(nK)
  
  for(i in 1:nK){
    d1 <- (log(S0/X[i])+(r-q+sig[i]*sig[i]/2)*capT)/sig[i]/sqrt(capT);
    d2 <- (log(S0/X[i])+(r-q-sig[i]*sig[i]/2)*capT)/sig[i]/sqrt(capT);

    pbar <- fh(d1,n)
    p <- fh(d2,n)
    u <- exp((r-q)*dt)*pbar/p
    d <- (exp((r-q)*dt)-p*u)/(1-p)
  
    # Loop over each node and calculate the CRR underlying price tree
    priceTree <- array(dim=c(steps+1,steps+1))
    priceTree[1,1] <- S0;
  
    for (idx in 2:(steps+1) )
    {
      priceTree[1:(idx-1),idx] <- priceTree[1:(idx-1),(idx-1)] * u;
      priceTree[idx,idx] <- priceTree[(idx-1),(idx-1)]*d;
    }
  
    # Calculate the value at expiry
    valueTree <- array(dim=c(steps+1,steps+1))
    
    if(oType == 'P')
    {
      valueTree[,] <- pmax(X[i]-priceTree[,],0)
    }
    else if(oType == 'C')
    {
      valueTree[,] <- pmax(priceTree[,]-X[i],0,na.rm = TRUE)
    }
    
    #Loop backwards to get values at the earlier times
    #steps = size(priceTree,2)-1;
    for (idx in steps:1)
    {
      valueTree[1:idx,idx] <- exp(-r*dt)*(p*valueTree[1:idx,(idx+1)] + (1-p)*valueTree[2:(idx+1),(idx+1)])
      if (earlyExercise == TRUE)
      {
        if(oType == 'P')
        {
          valueTree[1:idx,idx] <- pmax(X[i]-priceTree[1:idx,idx], valueTree[1:idx,idx], na.rm = TRUE)
        }
        else if(oType == 'C')
        {
          valueTree[1:idx,idx] <- pmax(priceTree[1:idx,idx]-X[i], valueTree[1:idx,idx], na.rm = TRUE)
        }
      }
    }

    V[i] <- (valueTree[1,1])
  }
  return (V)  
}


# This function now works with vectors of strikes and option values
binImpVolLR <- function(S0, K, capT, r, q,V, oType, earlyExercise, steps=10000)
{
  nK <- length(K);
  sigmaL <- rep(1e-10,nK);
  
  VL <- binPriceLR(K, S0, r, q,sigmaL, capT, oType, earlyExercise, steps)
  
  sigmaH <- rep(10,nK);
  VH <- binPriceLR(K, S0, r, q,sigmaH, capT, oType, earlyExercise, steps)
  
  while (mean(sigmaH - sigmaL) > 1e-10)
  {
    sigma <- (sigmaL + sigmaH)/2;
    VM <- binPriceLR(K, S0, r, q,sigma, capT, oType, earlyExercise, steps)
    
    VL <- VL + (VM < V)*(VM-VL);
    sigmaL <- sigmaL + (VM < V)*(sigma-sigmaL);
    VH <- VH + (VM >= V)*(VM-VH);
    sigmaH <- sigmaH + (VM >= V)*(sigma-sigmaH);
  }
  return(sigma);
}




