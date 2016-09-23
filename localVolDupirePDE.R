# Adapted from Rolf Poulsen's code
# In this version, sigma is a function of K and t
# Option prices of for all strikes and expirations returned from one call to the numerical PDE routine by solving the Dupire equation with appropriate boundary conditions.



fdLocalVolDupire <- function(S0,capT,r,sigma,dK,dt,sdwidth=4, theta=0.5, start=1)
{

delta <- 0; # We can reinstate dividendK later...
mu <- r-delta;

tridagsoln<-function(a,b,c,r)
{   n<-length(b)
    gam<-u<-tridagsoln<-1:n
    bet<-b[1]
    u[1]<-r[1]/bet
    for (j in 2:n) {
        gam[j]<-c[j-1]/bet
        bet<-b[j]-a[j]*gam[j]
        u[j]<-(r[j]-a[j]*u[j-1])/bet
    }
    for (j in (n-1):1) u[j]<-u[j]- gam[j+1]*u[j+1]

    tridagsoln<-u
}

tVec <- (1:10)/10*capT;
sigmaAvg <- mean(sigma(S0,tVec)); # Compute rough MLP estimate ATM

initialcond<-function(K){pmax(S0-K,0) };

highboundary<-function(K,texp){0};
lowboundary<-function(K,texp){S0*exp(-mu*texp)};

# set up finite-difference grid

Kmax<-S0*exp(sdwidth*sigmaAvg*sqrt(capT));
Kmin<-0;

nS0 <- ceiling(S0/dK);
dK <- S0/nS0; # Ensure that S0 is exactly on the grid

nK <- ceiling(Kmax/dK);
Kgrid<-(1:(nK-1))*dK;
Kmin <- 0;
Kmax <- nK*dK;

tvec<-capT
while (min(tvec) > 0){ tvec<-c(max(tvec[1]-dt,0),tvec)}

n.space<-length(Kgrid)
n.time<-length(tvec);

result<-matrix(nrow=(n.space+2),ncol=n.time); # Note that the whole grid is stored in memory!

KgridPlus <- c(Kmin,Kgrid,Kmax);
result[,n.time] <- initialcond(KgridPlus);
tm1 <- capT-dt+dt*theta; # Time computed consistently with implicitness
if (start==1) result[,(n.time-1)] <- BSFormula(S0,KgridPlus, dt, r, sigma(KgridPlus,tm1));

result[1,] <- lowboundary(Kmin,capT-tvec);
result[(n.space+2),] <- highboundary(Kmax,capT-tvec);

for (j in (n.time-1-start):1){

    dt <- tvec[j+1]-tvec[j];
    t1 <- tvec[j]+theta*dt; # Time chosen to be consistent with implicitness parameter theta
    
    # Note that these vectors are now time-dependent in general and need to be inside the time-loop
    vol <- sigma(Kgrid,t1);
    a <- ((1-theta)/(2*dK))*(mu*Kgrid-(vol*Kgrid)^2/dK);
    b <- 1/dt+(1-theta)*(r+(vol*Kgrid/dK)^2);
    c <- ((1-theta)/(2*dK))*(-mu*Kgrid-(vol*Kgrid)^2/dK);

    alpha <- ((-theta)/(2*dK))*(mu*Kgrid-(vol*Kgrid)^2/dK);
    beta <- 1/dt-theta*(r+(vol*Kgrid/dK)^2);
    gamma <- (theta/(2*dK))*(mu*Kgrid+(vol*Kgrid)^2/dK);
    eps <- 0;

    RHS <- alpha*result[1:n.space,(j+1)]+beta*result[2:(n.space+1),(j+1)]+gamma*result[3:(n.space+2),(j+1)]

    RHS[1] <- RHS[1]-a[1]*result[1,j]
    RHS[n.space] <- RHS[n.space]-c[n.space]*result[(n.space+2),j]

    result[2:(n.space+1),j] <- tridagsoln(a,b,c,RHS)
    
}

return(list(spacegrid=KgridPlus,timegrid=tvec,soln=result));

}


#res <- fdLocalVolDupire(S0=1,capT=1,r=0,sigma=function(x,t){.2},dK=0.001,dt=0.01,sdwidth=4, theta=0.5, start=1)

#------------------------------------------------------------------------------------------
# Generic version of PDE code
callLocalVolDupirePDE <- function(S0, r, q, sigma, t, dK, dt, sdw, start=1, theta=1/2){
    
        
    #tst <- fdLocalVol(S0,capT=t,strike=K,r,sigma,2*dK,dt,sdw,start=start,theta=theta);
    tst <- fdLocalVolDupire(S0,capT=t,r,sigma,2*dK,dt,sdw,start=start,theta=theta);
    Sout<-tst$spacegrid; 
    Sout <- Sout[abs(Sout-S0)<4*dK];
    w2h1k <- approx(tst$spacegrid,tst$soln[,1],Sout)$y

    #tst <- fdLocalVol(S0,capT=t,strike=K,r,sigma,dK,dt,sdw,start=start, theta=theta)
    tst <- fdLocalVolDupire(S0,capT=t,r,sigma,dK,dt,sdw,start=start, theta=theta)
    w1h1k <- approx(tst$spacegrid,tst$soln[,1],Sout)$y
    
    #tst <- fdLocalVol(S0,capT=t,strike=K,r,sigma,dK,2*dt,sdw,start=start,theta=theta)
    tst <- fdLocalVolDupire(S0,capT=t,r,sigma,dK,2*dt,sdw,start=start,theta=theta)
    w1h2k <- approx(tst$spacegrid,tst$soln[,1],Sout)$y

    wextrapol <- w1h1k+(1/3)*(w1h1k-w1h2k)+(1/3)*(w1h1k-w2h1k); # Richardson extrapolation
    return(approx(Sout,wextrapol,S0)$y);
    
    
    }
