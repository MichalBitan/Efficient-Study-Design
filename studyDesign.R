## @knitr sampleSizeNew

library("rgenoud")

#using optim for integers

#2 samples with integer optimization
#using the optim function finding optimum n and k for one group
sampleSize2groups<-function(sigmaEpsilon1, r_delta1, r_phi1, sigmaEpsilon2, r_delta2, r_phi2, C, CB, CQ)
{
  res<-genoud(optim2groupsI, nvars=5, data.type.int=TRUE, starting.values=c(floor(C*0.5/(CB+CQ)), 1, floor(C*0.5/(CB+CQ)), floor(C*0.5/(CB+CQ)), 1), Domains = matrix(c(4, 1, 4, 4, 1, floor(C/(CB+CQ)), 6, floor(C/(CB+CQ)), floor(C/(CB+CQ)), 6), nrow = 5, ncol = 2), additionalPar=c(sigmaEpsilon1, r_delta1, r_phi1, sigmaEpsilon2, r_delta2, r_phi2, C, CB, CQ), print.level=0, solution.tolerance=0)
  n1<-(res$par)[1]
  K1<-(res$par)[2]
  N1<-(res$par)[3]
  C1<-n1*K1*CB+N1*CQ
  ratioC<-round(C1/C, 2)
  ratios<-c(seq(ratioC-0.05, ratioC+0.05, by=0.01))
  g1<-NULL
  g2<-NULL
  SEs<-NULL
  for (i in 1:length(ratios))
  {
    C1<-ratios[i]*C
    C2<-C-ratios[i]*C
    g1<-rbind(g1, sampleSize1GroupI(C1, CB, CQ, sigmaEpsilon1, r_delta1, r_phi1, 4, 1))
    g2<-rbind(g2, sampleSize1GroupI(C2, CB, CQ, sigmaEpsilon2, r_delta2, r_phi2, 4, 1))
    SEs<-c(SEs, sqrt(g1[i,4]^2+g2[i,4]^2))
  }
  m<-which.min(SEs)
  
  sSize<-data.frame(n1=g1[m,1], N1=g1[m,2], K1=g1[m,3], n2=g2[m,1], N2=g2[m,2], K2=g2[m,3], ratioC=ratios[m], SE=SEs[m])
  
  return (sSize)
}

optim2groupsI<-function(par, additionalPar)
{
  sigmaEpsilon1<-additionalPar[1]
  r_delta1<-additionalPar[2]
  r_phi1<-additionalPar[3]
  sigmaEpsilon2<-additionalPar[4]
  r_delta2<-additionalPar[5]
  r_phi2<-additionalPar[6]
  C<-additionalPar[7]
  CB<-additionalPar[8]
  CQ<-additionalPar[9]
  n1<-par[1]
  K1<-par[2]
  N1<-par[3]
  n2<-par[4]
  K2<-par[5]
  
  N2<-floor(C/CQ-CB/CQ*(n1*K1+n2*K2)-N1)
  a1<-(K1+r_delta1)*(N1*n1-2*N1-n1)/K1
  b1<-(n1-2)*(N1-n1)/(1+r_phi1)
  a2<-(K2+r_delta2)*(N2*n2-2*N2-n2)/K2
  b2<-(n2-2)*(N2-n2)/(1+r_phi2)
  
  if (n1*K1*CB+N1*CQ+n2*K2*CB+N2*CQ>C | N1<n1 | a1<b1 | N2<n2 | a2<b2)
    res<-100000000000
  else
    res<-sqrt(MuRatioSE(sigmaEpsilon1, r_delta1, r_phi1, n1, K1, N1)^2+MuRatioSE(sigmaEpsilon2, r_delta2, r_phi2, n2, K2, N2)^2)
  
  return (res)
}


#calculating the samples design study when using the error alpha the power and the diff between the two groups
#2 samples with integer optimization
#using the optim function findig optimum n and k for one group
studyDesign2groups<-function(error1, power, delta, sigmaEpsilon1, r_delta1, r_phi1, sigmaEpsilon2, r_delta2, r_phi2, CB, CQ)
{
  SE<-delta/(qnorm(1-error1/2)+qnorm(power))
  
  r_K1<-r_delta1/(CB/CQ)
  r_K2<-r_delta2/(CB/CQ)
  if(r_K1<2)
    K1<-1
  else if (r_K1<6 & r_K1>=2)
    K1<-2
  else if (r_K1<12 & r_K1>=6)
    K1<-3
  else
    K1<-4
  if(r_K2<2)
    K2<-1
  else if (r_K2<6 & r_K2>=2)
    K2<-2
  else if (r_K2<12 & r_K2>=6)
    K2<-3
  else
    K2<-4
  
  C<-2/SE^2*((sigmaEpsilon1*(1+r_delta1/K1)*(CQ+K1*CB))+(sigmaEpsilon2*(1+r_delta2/K2)*(CQ+K2*CB)))

  N1<-round(0.5*C/(CQ+K1*CB), 0)
  N2<-round(0.5*C/(CQ+K2*CB), 0)
  SE_C0<-sqrt(MuRatioSE(sigmaEpsilon1, r_delta1, r_phi1, N1, K1, N1)^2+MuRatioSE(sigmaEpsilon2, r_delta2, r_phi2, N2, K2, N2)^2)
  res<-data.frame(n1=N1, N1=N1, K1=K1, n2=N2, N2=N2, K2=K2, ratioC=0.5, SE=SE_C0)
  res<-sampleSize2groups(sigmaEpsilon1, r_delta1, r_phi1, sigmaEpsilon2, r_delta2, r_phi2, C, CB, CQ)
  
  for(i in 1:100)
  {
    C<-C*(res$SE/SE)^2
    res<-sampleSize2groups(sigmaEpsilon1, r_delta1, r_phi1, sigmaEpsilon2, r_delta2, r_phi2, C, CB, CQ)
    if (res$SE<=SE+tol & res$SE >=SE-tol)
    {
      res<-cbind(res, C=C)
      return (res)
    }
  }
  
}


#with tol=tolerence how much far from the target SE
studyDesign2groups2<-function(error1, power, delta, sigmaEpsilon1, r_delta1, r_phi1, sigmaEpsilon2, r_delta2, r_phi2, CB, CQ, tol=0.001)
{
  SE<-delta/(qnorm(1-error1/2)+qnorm(power))
  
  r_K1<-r_delta1/(CB/CQ)
  r_K2<-r_delta2/(CB/CQ)
  if(r_K1<2)
    K1<-1
  else if (r_K1<6 & r_K1>=2)
    K1<-2
  else if (r_K1<12 & r_K1>=6)
    K1<-3
  else
    K1<-4
  if(r_K2<2)
    K2<-1
  else if (r_K2<6 & r_K2>=2)
    K2<-2
  else if (r_K2<12 & r_K2>=6)
    K2<-3
  else
    K2<-4
  
  
  C<-2/SE*sqrt((sigmaEpsilon1*(CQ*r_delta1/K1+K1*CB+CQ+CB*r_delta1))^2+(sigmaEpsilon2*(CQ*r_delta1/K2+K2*CB+CQ+CB*r_delta2))^2)
  res<-sampleSize2groups(sigmaEpsilon1, r_delta1, r_phi1, sigmaEpsilon2, r_delta2, r_phi2, C, CB, CQ)
  if (res$SE<=SE+tol & res$SE >=SE-tol)
    return (res)
  
  for(i in 1:100)
  {
    C<-C*(res$SE/SE)^2
    res<-sampleSize2groups(sigmaEpsilon1, r_delta1, r_phi1, sigmaEpsilon2, r_delta2, r_phi2, C, CB, CQ)
    if (res$SE<=SE+tol & res$SE >=SE-tol)
    {
      res<-cbind(res, C=C)
      return (res)
    }
  }
  
}



sampleSize2groupsSameK<-function(sigmaEpsilon1, r_delta1, r_phi1, sigmaEpsilon2, r_delta2, r_phi2, C, CB, CQ)
{
  res<-genoud(optim2groupsISameK, nvars=4, data.type.int=TRUE, starting.values=c(floor(C*0.5/(CB+CQ)), 1, floor(C*0.5/(CB+CQ)), floor(C*0.5/(CB+CQ))), Domains = matrix(c(4, 1, 4, 4, floor(C/(CB+CQ)), 6, floor(C/(CB+CQ)), floor(C/(CB+CQ))), nrow = 4, ncol = 2), additionalPar=c(sigmaEpsilon1, r_delta1, r_phi1, sigmaEpsilon2, r_delta2, r_phi2, C, CB, CQ), print.level=0, solution.tolerance=0)
  n1<-(res$par)[1]
  K<-(res$par)[2]
  N1<-(res$par)[3]
  n2<-(res$par)[4]
  C1<-n1*K*CB+N1*CQ
  ratioC<-round(C1/C, 2)
  N2<-floor(C/CQ-CB/CQ*K*(n1+n2)-N1)

    sSize<-data.frame(n1=n1, N1=N1, K1=K, n2=n2, N2=N2, K2=K, ratioC=ratioC, SE=res$value)
  
  return (sSize)
}




#the optimum function to find n with fixed K
optim2groupsISameK<-function(par, additionalPar)
{
  sigmaEpsilon1<-additionalPar[1]
  r_delta1<-additionalPar[2]
  r_phi1<-additionalPar[3]
  sigmaEpsilon2<-additionalPar[4]
  r_delta2<-additionalPar[5]
  r_phi2<-additionalPar[6]
  C<-additionalPar[7]
  CB<-additionalPar[8]
  CQ<-additionalPar[9]
  n1<-par[1]
  K1<-par[2]
  N1<-par[3]
  n2<-par[4]
  K2<-par[2]
  
  N2<-floor(C/CQ-CB/CQ*(n1*K1+n2*K2)-N1)
  a1<-(K1+r_delta1)*(N1*n1-2*N1-n1)/K1
  b1<-(n1-2)*(N1-n1)/(1+r_phi1)
  a2<-(K2+r_delta2)*(N2*n2-2*N2-n2)/K2
  b2<-(n2-2)*(N2-n2)/(1+r_phi2)
  
  if (n1*K1*CB+N1*CQ+n2*K2*CB+N2*CQ>C | N1<n1 | a1<b1 | N2<n2 | a2<b2)
    res<-100000000000
  else
    res<-sqrt(MuRatioSE(sigmaEpsilon1, r_delta1, r_phi1, n1, K1, N1)^2+MuRatioSE(sigmaEpsilon2, r_delta2, r_phi2, n2, K2, N2)^2)
  
  return (res)
}






MuRatioSE<-function(sigmaEpsilon, r_delta, r_phi, n, K, N)
{
  a<-(K+r_delta)*(N*n-2*N-n)/K
  b<-(n-2)*(N-n)/(1+r_phi)
  SE<-sqrt(sigmaEpsilon/(N*n*(n-3))*(a-b))
  return (SE)
}



#using the optim function findig optimum n and k for one group
sampleSize1group<-function(C, CB, CQ, sigmaEpsilon, r_delta, r_phi, minn, minK)
{
  sSize<-NULL
  N<-(C-minn*minK*CB)/CQ
  if (minn<4)
    minn<-4
  for (K in minK:6)
  {
    if (C/(CQ+K*CB)<minn)
      break
    res<-optimize(f=optimn,lower=minn, upper=C/(CQ+K*CB), sigmaEpsilon=sigmaEpsilon, r_delta=r_delta, r_phi=r_phi, C=C, CB=CB, CQ=CQ, K=K)
    n<-res$minimum
    N<-(C-n*K*CB)/CQ
    sSize<-rbind(sSize, data.frame(n=n, N=N, K=K, SE=res$objective))
  }
  minSE<-sSize[which.min(sSize$SE),]
  
  return (minSE)
}




#one sample with integer optimization
#using the optim function finding optimum n and k for one group
sampleSize1GroupI<-function(C, CB, CQ, sigmaEpsilon, r_delta, r_phi, minn, minK)
{
  sSize<-NULL
  N<-floor((C-minn*minK*CB)/CQ)
  if (minn<4)
    minn<-4
  res<-genoud(optimn, nvars=2, data.type.int=TRUE, starting.values=c(C/(CQ+CB*minK)*0.7, minK), Domains = matrix(c(minn, minK, floor(C/(CQ+minK*CB)), 6), nrow = 2, ncol = 2), sigmaEpsilon=sigmaEpsilon, r_delta=r_delta, r_phi=r_phi, C=C, CB=CB, CQ=CQ, print.level=0)
    n<-(res$par)[1]
    K<-(res$par)[2]
    N<-floor((C-n*K*CB)/CQ)
    sSize<-rbind(sSize, data.frame(n=n, N=N, K=K, SE=res$value))

  return (sSize)
}

#the optimum function to find n with fixed K
optimn<-function(par, sigmaEpsilon, r_delta, r_phi, C, CB, CQ)
{
  n<-par[1]
  K<-par[2]
  N<-floor((C-n*K*CB)/CQ)
  a<-(K+r_delta)*(N*n-2*N-n)/K
  b<-(n-2)*(N-n)/(1+r_phi)
  
  if (n*K*CB+N*CQ>C | N<n | a<b)
    res<-100000000000
  else
    res<-MuRatioSE(sigmaEpsilon, r_delta, r_phi, n, K, N)
  
  return (res)
}


#one sample with integer optimization
#using the optim function finding optimum n for one group when K is fixed
sampleSize1GroupIfixedK<-function(C, CB, CQ, sigmaEpsilon, r_delta, r_phi, minn, K)
{
  sSize<-NULL
  N<-floor((C-minn*K*CB)/CQ)
  if (minn<4)
    minn<-4
  res<-genoud(optimnFixedK, nvars=1, data.type.int=TRUE, starting.values=c(C/(CQ+CB*K)*0.7), Domains = matrix(c(minn, floor(C/(CQ+K*CB))), nrow = 1, ncol = 2), sigmaEpsilon=sigmaEpsilon, r_delta=r_delta, r_phi=r_phi, C=C, CB=CB, CQ=CQ, K=K, print.level=0)
  n<-(res$par)[1]
  N<-floor((C-n*K*CB)/CQ)
  sSize<-rbind(sSize, data.frame(n=n, N=N, K=K, SE=res$value))
  
  return (sSize)
}

#the optimum function to find n with fixed K
optimnFixedK<-function(par, sigmaEpsilon, r_delta, r_phi, C, CB, CQ, K)
{
  n<-par[1]
  N<-floor((C-n*K*CB)/CQ)
  a<-(K+r_delta)*(N*n-2*N-n)/K
  b<-(n-2)*(N-n)/(1+r_phi)
  
  if (n*K*CB+N*CQ>C | N<n | a<b)
    res<-100000000000
  else
    res<-MuRatioSE(sigmaEpsilon, r_delta, r_phi, n, K, N)
  
  return (res)
}



#when n=N find optim K with integer
#one sample with integer optimization
sampleSize1GroupInN<-function(C, CB, CQ, sigmaEpsilon, r_delta)
{
  sSize<-NULL
  res<-genoud(optimnnN, nvars=1, data.type.int=TRUE, starting.values=c(1), Domains = matrix(c(1, C/CB), nrow = 1, ncol = 2), sigmaEpsilon=sigmaEpsilon, r_delta=r_delta, C=C, CB=CB, CQ=CQ, print.level=0)
  n<-C/(CQ+K*CB)
  N<-n
  sSize<-rbind(sSize, data.frame(n=n, N=N, K=K, SE=res$value))
  
  return (sSize)
}

#the optimum function to find K with n=N
optimnnN<-function(par, sigmaEpsilon, r_delta, r_phi, C, CB, CQ, K)
{
  K<-par[1]
  N<-C/(CQ+K*CB)

  res<-sigmaEpsilon*(CQ+K*CB)/C*(1+r_delta/K)
  
  return (res)
}


sampleSizenN<-function(C, CB, CQ, sigmaEpsilon, r_delta)
{
  K<-CQ*r_delta/(sqrt(CB))
  N<-C/(CQ+K*CB)
  SE<-sigmaEpsilon*(CQ+K*CB)/C*(1+r_delta/K)
  ssize<-data.frame(n=N, N=N, K=K, SE=SE)
  return (ssize)
}

