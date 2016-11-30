# Adapted from Rolf Poulsen's code
# In this version, sigma is a function of S and t

fdLocalVol <- function(S0,capT,strike,r,sigma,dS,dt,sdwidth=4, theta=0.5, start=1,oType="C")
{
delta <- 0; # We can reinstate dividends later...
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

#sVec <- S0+(strike-S0)*(1:10)/10;
#tVec <- (1:10)/10*capT;

### need to change the following line (local vol function cannot take t=0)
#sigmaAvg <- (sigma(S0,0)+sigma(strike,T))/2; # Compute rough MLP estimate
sigmaAvg <- sigma(log(S0/S0),capT); # Compute rough MLP estimate

if(oType=="C"){
  initialcond<-function(S){pmax(S-strike,0) };
  highboundary<-function(S,timeleft){pmax(S*exp(-delta*timeleft)-exp(-r*timeleft)*strike,S-strike)};
  lowboundary<-function(S,timeleft){0};
}
else{
  initialcond<-function(S){pmax(strike-S,0) };
  highboundary<-function(S,timeleft){0};
  lowboundary<-function(S,timeleft){pmax(exp(-r*timeleft)*strike,strike-S)};
}


# set up finite-difference grid

Smax<-S0*(1+sdwidth*sigmaAvg*sqrt(capT));
Smin<-max(2*dS,S0*(1-sdwidth*sigmaAvg*sqrt(capT)));

Sgrid<-strike;

while (max(Sgrid) < Smax){Sgrid<-c(Sgrid,max(Sgrid)+dS)}
while (min(Sgrid) > Smin){Sgrid<-c(min(Sgrid)-dS,Sgrid)}

tvec<-capT
while (min(tvec) > 0){ tvec<-c(max(tvec[1]-dt,0),tvec)}

Smin<-min(Sgrid)-dS
Smax<-max(Sgrid)+dS

n.space<-length(Sgrid);
n.time<-length(tvec);

result<-matrix(nrow=(n.space+2),ncol=n.time); # Note that the whole grid is stored in memory!

SGridPlus <- c(Smin,Sgrid,Smax);
result[,n.time] <- initialcond(SGridPlus);
tm1 <- capT-dt+dt*theta; # Time computed consistently with implicitness
if (start==1){
  if(oType=="C"){
    result[,(n.time-1)] <- pmax(BSFormula(SGridPlus, strike, dt, r, sigma(log(strike/SGridPlus),capT-tm1)),SGridPlus-strike);
  }
  else{
    result[,(n.time-1)] <- pmax(BSFormulaPut(SGridPlus, strike, dt, r, sigma(log(strike/SGridPlus),capT-tm1)),strike-SGridPlus);
  }
}
  
result[1,] <- lowboundary(Smin,capT-tvec);
result[(n.space+2),] <- highboundary(Smax,capT-tvec);

for (j in (n.time-1-start):1){

    dt <- tvec[j+1]-tvec[j];
    t1 <- tvec[j]+theta*dt; # Time chosen to be consistent with implicitness parameter theta
    
    # Note that these vectors are now time-dependent in general and need to be inside the time-loop
    vol <- sigma(log(strike/Sgrid),capT-t1);
    a <- ((1-theta)/(2*dS))*(mu*Sgrid-(vol*Sgrid)^2/dS);
    b <- 1/dt+(1-theta)*(r+(vol*Sgrid/dS)^2);
    c <- ((1-theta)/(2*dS))*(-mu*Sgrid-(vol*Sgrid)^2/dS);

    alpha <- ((-theta)/(2*dS))*(mu*Sgrid-(vol*Sgrid)^2/dS);
    beta <- 1/dt-theta*(r+(vol*Sgrid/dS)^2);
    gamma <- (theta/(2*dS))*(mu*Sgrid+(vol*Sgrid)^2/dS);
    eps <- 0;

    RHS <- alpha*result[1:n.space,(j+1)]+beta*result[2:(n.space+1),(j+1)]+gamma*result[3:(n.space+2),(j+1)]

    RHS[1] <- RHS[1]-a[1]*result[1,j]
    RHS[n.space] <- RHS[n.space]-c[n.space]*result[(n.space+2),j]

    result[2:(n.space+1),j] <- tridagsoln(a,b,c,RHS)
    if(oType=="C"){
      result[2:(n.space+1),j]=pmax(result[2:(n.space+1),j],Sgrid-strike)
    }
    else{
      result[2:(n.space+1),j]=pmax(result[2:(n.space+1),j],strike-Sgrid)
    }
}

return(list(spacegrid=SGridPlus,timegrid=tvec,soln=result));

}

#------------------------------------------------------------------------------------------
# Generic version of PDE code
LocalVolPDEAmerican <- function(S0, K, r, q, sigma, t, dS, dt, sdw, start=1, theta=1/2,oType="C"){
    
        
    tst <- fdLocalVol(S0,capT=t,strike=K,r,sigma,2*dS,dt,sdw,start=start,theta=theta,oType);
    Sout<-tst$spacegrid; 
    Sout <- Sout[abs(Sout-S0)<4*dS];
    w2h1k <- approx(tst$spacegrid,tst$soln[,1],Sout)$y

    tst <- fdLocalVol(S0,capT=t,strike=K,r,sigma,dS,dt,sdw,start=start, theta=theta,oType)
    w1h1k <- approx(tst$spacegrid,tst$soln[,1],Sout)$y
    
    tst <- fdLocalVol(S0,capT=t,strike=K,r,sigma,dS,2*dt,sdw,start=start,theta=theta,oType)
    w1h2k <- approx(tst$spacegrid,tst$soln[,1],Sout)$y

    wextrapol <- w1h1k+(1/3)*(w1h1k-w1h2k)+(1/3)*(w1h1k-w2h1k); # Richardson extrapolation
    
    return(approx(Sout,wextrapol,S0)$y);
    
    
    }
