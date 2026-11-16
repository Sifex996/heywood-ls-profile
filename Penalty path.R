# Penalty path.R
# Sect. 6.1 and 6.7: the penalty path of the 3 x 3 boundary example (Table 1), and the numerical
# checks of the results on the family P_w(L) = Q(L) + w/2 sum_i (C_ii - S_ii)_+^2.

source("Estimators.R")

set.seed(123)

# Table 1 ------------------------------------------------------------------------------------
S=matrix(c(1,0.60,0.60,0.60,1,0.25,0.60,0.25,1),3,3)
t=uniroot(function(t) t^3+0.75*t-0.60,c(0,1),tol=1e-15)$root
M_H=2*(t-0.60)^2+(t^2-0.25)^2
mu_bar=2*t*(0.60-t)

# P_w on a grid of lambda in [0, 1.6]^3, refined around the best point
Pw_grid=function(w)
{
  lo=rep(0,3)
  hi=rep(1.6,3)
  for (it in 1:40)
  {
    g=expand.grid(seq(lo[1],hi[1],length.out=41),seq(lo[2],hi[2],length.out=41),seq(lo[3],hi[3],length.out=41))
    Q=(0.60-g[,1]*g[,2])^2+(0.60-g[,1]*g[,3])^2+(0.25-g[,2]*g[,3])^2
    P=Q+0.5*w*(pmax(g[,1]^2-1,0)^2+pmax(g[,2]^2-1,0)^2+pmax(g[,3]^2-1,0)^2)
    best=unlist(g[which.min(P),])
    h=(hi-lo)/8
    lo=pmax(best-h,0)
    hi=pmin(best+h,1.6)
  }
  lamda=best
  Q=(0.60-lamda[1]*lamda[2])^2+(0.60-lamda[1]*lamda[3])^2+(0.25-lamda[2]*lamda[3])^2
  return(c(M=min(P),Q=Q,r1=max(lamda[1]^2-1,0),lamda))
}

wval=c(0.1,0.5,1,2,10,100)
Tab1=matrix(NA,6,6,dimnames=list(NULL,c("w","M(w)","Q_w","r_1","bound","w r_1")))
for (k in 1:6)
{
  x=Pw_grid(wval[k])
  Tab1[k,]=c(wval[k],x[1],x[2],x[3],sqrt(2*(M_H-x[2])/wval[k]),wval[k]*x[3])
}
round(Tab1,6)
c(t=t,M_H=M_H,mu_bar=mu_bar)

# Checks of Sect. 6.7 on 322 matrices -----------------------------------------------------------
wgrid=c(0.05,0.2,0.5,1,2,5,20,100,1000)
tol=1e-7

# Best of several BFGS starts on P_w
best_Pw=function(S,w,starts,rtol=reltol)
{
  best=NULL
  for (j in 1:length(starts))
  {
    f=tryCatch(fit_Pw(S,starts[[j]],w,rtol),error=function(e) NULL)
    if (!is.null(f) && (is.null(best) || f$value<best$value)) best=f
  }
  best$C=tcrossprod(best$Lamda)
  best$Q=Qval(S,best$C)
  best$v=Pen(S,best$C)
  best$r=pmax(diag(best$C)-diag(S),0)
  return(best)
}

check=function(S,m,nst)
{
  e1=fit_E1(S,m,nst)
  e2=fit_E2(S,m,nst)
  L1=lamda_C(e1$C,m)
  L2=lamda_C(e2$C,m)
  M0=e1$Q
  MH=e2$Q
  v0=Pen(S,e1$C)
  st1=sqrt(sum((-2*off(S-tcrossprod(L1))%*%L1)^2))<1e-6
  g2=sqrt(sum((-2*off(S-tcrossprod(L2))%*%L2)^2))
  starts=list()
  for (j in 1:nst)
  {
    starts[[j]]=start_L(S,m,j)
  }
  starts=c(starts,list(L1,L2))
  # Forward sweep in w with warm starts, then a backward sweep
  fits=list()
  for (k in 1:length(wgrid))
  {
    if (k==1) st=starts else st=c(starts,list(fits[[k-1]]$Lamda))
    fits[[k]]=best_Pw(S,wgrid[k],st)
  }
  for (k in (length(wgrid)-1):1)
  {
    alt=best_Pw(S,wgrid[k],list(fits[[k+1]]$Lamda))
    if (alt$value<fits[[k]]$value) fits[[k]]=alt
  }
  Qw=sapply(fits,function(f) f$Q)
  vw=sapply(fits,function(f) f$v)
  Mw=sapply(fits,function(f) f$value)
  rmax=sapply(fits,function(f) max(f$r))
  rowb=sqrt(2*(rowSums(S^2)-diag(S)^2))
  A5=TRUE
  for (k in 1:length(wgrid))
  {
    A5=A5 && all(fits[[k]]$r<=rowb/sqrt(wgrid[k])+tol)
  }
  interior=MH-M0<=tol
  binding=g2>1e-5
  return(c(A1_Q=all(diff(Qw)>=-tol),A1_v=all(diff(vw)<=tol),A2=all(Mw>=M0-tol & Mw<=MH+tol),
           A3=all(vw<=(MH-Qw)/wgrid+tol),A4=ifelse(st1,all(vw<=v0+tol),NA),A5=A5,
           A6_interior=ifelse(interior,all(rmax<=1e-5),NA),
           A6_binding=ifelse(binding,all(rmax>0 & Mw>M0+1e-12 & Mw<MH-1e-12),NA),
           interior=interior,binding=binding,gap=(MH-Qw[9])/max(MH,1e-12),r1=rmax[4],
           bound=sqrt(2*max(MH-Qw[4],0)),bound_old=sqrt(2*sum(S[upper.tri(S)]^2))))
}

# The 46 populations of the simulation, and three samples from each at n = 50 and n = 200
keys=c("block","m","k","pattern","truth","lambda1")
A=expand.grid(truth=c("T0","T1","M1","M2","M2s"),pattern=c("moderate","high"),k=c(3,5),m=c(1,3),
              stringsAsFactors=FALSE)
A$block="A"
A$lambda1=NA
B=expand.grid(lambda1=c(1.1,1.2,1.3),k=c(3,6),stringsAsFactors=FALSE)
B$m=1
B$truth="B"
B$pattern="overshoot"
B$block="B"
pops=rbind(A[,keys],B[,keys])

Chk=NULL
for (q in 1:nrow(pops))
{
  pop=build_pop(pops$block[q],pops$m[q],pops$k[q],pops$pattern[q],pops$truth[q],pops$lambda1[q])
  m=pops$m[q]
  if (m==3) nst=10 else nst=5
  Chk=rbind(Chk,check(pop$Sigma,m,nst))
  for (n in c(50,200))
  {
    for (rep in 1:3)
    {
      Chk=rbind(Chk,check(draw_S(chol(pop$Sigma),n),m,nst))
    }
  }
}
nrow(Chk)
# Applicable cases, and cases in which the statement holds
rbind(applicable=colSums(!is.na(Chk[,1:8])),holds=colSums(Chk[,1:8]==1,na.rm=TRUE))
c(sum(Chk[,"interior"]),sum(Chk[,"binding"]))
max(Chk[Chk[,"binding"]==1,"gap"])
pos=Chk[,"r1"]>1e-6
c(sum(pos),median(Chk[pos,"r1"]),median(Chk[pos,"bound"]),median(Chk[pos,"bound_old"]))
range(Chk[pos,"bound"]/Chk[pos,"r1"])

# Multiplier of the cap (Proposition on the multiplier) ------------------------------------------
# mu_i = -(grad Q(L)_i . l_i) / (2 |l_i|^2) on the active rows of the capped fit
wmu=c(1,10,100,1000)
Mu=NULL
for (q in 1:nrow(pops))
{
  pop=build_pop(pops$block[q],pops$m[q],pops$k[q],pops$pattern[q],pops$truth[q],pops$lambda1[q])
  m=pops$m[q]
  if (m==3) nst=10 else nst=5
  for (s in 1:2)
  {
    if (s==1) S=pop$Sigma else S=draw_S(chol(pop$Sigma),200)
    e2=fit_E2(S,m,nst)
    L2=lamda_C(e2$C,m)
    G=-2*off(S-tcrossprod(L2))%*%L2
    mu=-rowSums(G*L2)/(2*rowSums(L2^2))
    active=diag(e2$C)>=diag(S)*(1-1e-6)
    if (any(active))
    {
      i=which(active)[which.max(mu[active])]
      starts=list()
      for (j in 1:nst)
      {
        starts[[j]]=start_L(S,m,j)
      }
      starts=c(starts,list(L2))
      wr=vector()
      prev=NULL
      for (k in 1:4)
      {
        f=best_Pw(S,wmu[k],c(starts,list(prev)),1e-15)
        wr[k]=wmu[k]*f$r[i]
        prev=f$Lamda
      }
      Mu=rbind(Mu,c(mu=mu[i],wr))
    }
  }
}
Mu=Mu[Mu[,1]>1e-6,]
nrow(Mu)
max(abs(Mu[,5]-Mu[,1])/Mu[,1])
sum(Mu[,2]<=Mu[,3]+1e-8 & Mu[,3]<=Mu[,4]+1e-8 & Mu[,4]<=Mu[,5]+1e-6)
c(median(abs(Mu[,2]-Mu[,1])/Mu[,1]),max(abs(Mu[,2]-Mu[,1])/Mu[,1]))
