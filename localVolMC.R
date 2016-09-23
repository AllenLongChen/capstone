evolveEuler <- function(x,capT,dt,Z,i){
  
  #local variance
  t<-capT-(i-1)*dt
  v<-a+b*(rho*(-x/sqrt(t)-m1)+sqrt((-x/sqrt(t)-m1)^2+sig^2*t))
  
  # Log-stock process
  x <- x  - v/2*dt + sqrt(v*dt) * Z
  # Impose martingale constraint
  x <- x - log(mean(exp(x))) 
  
  return(x)
}

is.even <- function(j){as.logical((j+1) %% 2)} # A little function needed later

LocalVolMC <- function(params){
  
  res <- function(S0, capT, AK, N, m, evolve,exactVols=NULL)
  {
    
    a<<-params$a
    b<<-params$b
    sig<<-params$sig
    rho<<-params$rho
    m1<<-params$m1
    
    n <- m*2  #n is number of timesteps = 2*m so we can use Richardson extrapolation
    sqrt2 <- sqrt(2) 
    
    # We use a vertical array, one element per M.C. path
    x <- rep(0,N);
    xm <- x;
    W1m <- rep(0,N);
    
    # Loop for bias computation (N small, n big)
    for (i in 1:n)
    {
      # Two sets of correlated normal random vars.
      
      W1 <- rnorm(N) 
      W1 <- W1 - mean(W1);  W1 <- W1/sd(W1) 
      # Now W1 and W2 are forced to have mean=0 and sd=1
      
      # Add code for subgrid
      W1m <- W1m + W1/sqrt2;# N(0,1) rv's for subgrid
      
      if (is.even(i)) {
        #print(c(i,mean(W1m),mean(W2m),sd(W1m),sd(W2m),cor(W1m,W2m))) 
        resm <- evolve(xm,capT,capT/m,W1m,i/2) 
        xm <- resm
        W1m <- rep(0,N);
      }
      
      res <- evolve(x,capT,capT/n,W1,i) 
      x <- res
      
    }
    S <- S0*exp(x) 
    Sm <- S0*exp(xm) 
    
    # Now we have M vectors of final stock prices
    
    M <- length(AK) 
    AV <- numeric(M); AVdev <- numeric(M) 
    BSV <- numeric(M); BSVH <- numeric(M);  BSVL <- numeric(M) 
    iv2SD <- numeric(M);  bias <- numeric(M) 
    AVm <- numeric(M);  AVmdev <- numeric(M) 
    BSVm <- numeric(M);  BSVHm <- numeric(M);  BSVLm <- numeric(M) 
    iv2SDm <- numeric(M) 
    AV1 <- numeric(M);  AV1dev <- numeric(M) 
    BSV1 <- numeric(M); BSVH1 <- numeric(M); BSVL1 <- numeric(M) 
    iv2SDrom <- numeric(M); biasRom <- numeric(M)
    
    # Evaluate mean call value for each path
    for (i in 1:M)
    {
      # 2*m timesteps
      K <- AK[i] 
      V <- (S>K)*(S - K)  # Boundary condition for European call
      AV[i] <- mean(V)
      AVdev[i] <- sqrt(var(V)/length(V)) 
      
      BSV[i] <- BSImpliedVolCall(S0, K, capT, 0, AV[i]) 
      BSVL[i] <- BSImpliedVolCall(S0, K, capT, 0, AV[i] - AVdev[i]) 
      BSVH[i] <- BSImpliedVolCall(S0, K, capT, 0, AV[i] + AVdev[i]) 
      iv2SD[i] <- (BSVH[i]-BSVL[i]) 
      
      # m timesteps
      Vm <- (Sm>K)*(Sm - K)  # Boundary condition for European call
      AVm[i] <- mean(Vm) 
      AVmdev[i] <- sd(Vm) / sqrt(N) 
      BSVm[i] <- BSImpliedVolCall(S0, K, capT, 0, AVm[i]) 
      BSVLm[i] <- BSImpliedVolCall(S0, K, capT, 0, AVm[i] - AVmdev[i]) 
      BSVHm[i] <- BSImpliedVolCall(S0, K, capT, 0, AVm[i] + AVmdev[i]) 
      iv2SDm[i] <- (BSVH[i]-BSVL[i]) 
      
      # Romberg estimates 
      V1 <- 2*V - Vm 
      AV1[i] <- mean(V1) 
      AV1dev[i] <- sd(V1) / sqrt(N) 
      BSV1[i] <- BSImpliedVolCall(S0, K, capT, 0, AV1[i]) 
      BSVL1[i] <- BSImpliedVolCall(S0, K, capT, 0, AV1[i] - AV1dev[i]) 
      BSVH1[i] <- BSImpliedVolCall(S0, K, capT, 0, AV1[i] + AV1dev[i]) 
      iv2SDrom[i] <- (BSVH1[i]-BSVL1[i]) 
      
      if(!is.null(exactVols)) {bias <- BSV-exactVols} 
      if(!is.null(exactVols)) {biasRom <- BSV1-exactVols} 
    }
    
    l.AK <- length(AK)      
    data.out <- data.frame(AK,rep(N,l.AK),rep(2*m,l.AK),BSV,bias,iv2SD,BSVm,BSV1,biasRom,iv2SDrom) 
    names(data.out) <- c("Strikes","Paths","Steps","ivol","bias","twoSd","ivolm", "ivolRichardson", "biasRichardson","twoSdRichardson") 
    return(data.out) 
    
    
  }
  return(res) 
}