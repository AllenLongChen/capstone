source("BlackScholes.R");
source("localVolPDE.R");

# Example Andersen Quadratic model
lvA <- function(s,t){
    
    sig <- 0.2; # Rough parameters from Andersen (2008)
    psi <- -0.5;
    gamma <- .05;
    s0 <- 1;
    vol <- sig/s*(psi*s+(1-psi*s0)+gamma/2*(s-s0)^2/s0);
    return(vol);

}

# Compute option values using PDE
callValue <- function(k){sapply(k,function(k){callLocalVolPDE(S0=1, K=k, r=0, q=0, sigma=lvA, t=1, dS=0.01, dt=0.01, sdw=15)})};

# Graph PDE solution
k <- seq(0.5,1.5,0.05);
impvol <- function(k){BSImpliedVolCall(1,k,1,0,callValue(k))};
volk <- impvol(k);
plot(k,volk,col="red",type="b")
