# Row-wise algorithm benchmark.R
# Online Resource 1, Sect. S3: profiled block-coordinate descent (PBCD) for ULS with nonnegative
# uniquenesses against L-BFGS-B over psi (two full eigendecompositions per point, one, or a partial
# Lanczos eigensolver), BFGS over the loadings, and psych::fa(fm = "uls") for timing context.
# Every method starts from the same five starting values per matrix.

library("psych")
library("RSpectra")
source("Estimators.R")

set.seed(123)
nst=5

# P(L) and its gradient
pbcd_value=function(S,Lamda)
{
  C=tcrossprod(Lamda)
  return(Qval(S,C)+Pen(S,C))
}

pbcd_grad=function(S,Lamda)
{
  C=tcrossprod(Lamda)
  e=pmax(diag(C)-diag(S),0)
  return(-2*off(S-C)%*%Lamda+2*e*Lamda)
}

# Root of phi(mu) = sum b^2/(d+mu)^2 - s - mu on (0, Inf): safeguarded Newton
secular=function(d,b,s)
{
  phi=function(mu)
  {
    return(sum(b^2/(d+mu)^2)-s-mu)
  }
  lo=0
  hi=max(1,sqrt(sum(b^2))^(2/3))
  while (phi(hi)>0) hi=2*hi
  if (all(d>0)) mu=0 else mu=min(hi,1e-12*hi)
  for (it in 1:200)
  {
    f=phi(mu)
    if (abs(f)<=1e-15*max(1,s)) break
    if (f>0) lo=mu else hi=mu
    step=mu-f/(-2*sum(b^2/(d+mu)^3)-1)
    if (is.finite(step) && step>lo && step<hi) mu=step else mu=0.5*(lo+hi)
    if (hi-lo<=1e-15*max(1,hi)) break
  }
  return(mu)
}

# Exact row update: min ||s - A x||^2 + 1/2 (||x||^2 - s_ii)_+^2 (+ tau/2 ||x - x_old||^2)
row_update=function(H,g,s,tau,xold)
{
  m=length(g)
  if (tau>0)
  {
    H=H+diag(tau,m)
    g=g+tau*xold
  }
  if (all(g==0)) return(c(numeric(m),0))
  eg=eigen(H,symmetric=TRUE)
  d=pmax(eg$values,0)
  U=eg$vectors
  b=drop(crossprod(U,g))
  pos=d>max(d)*1e-12
  b[!pos]=0
  c0=numeric(m)
  c0[pos]=b[pos]/d[pos]
  if (sum(c0^2)<=s) return(c(drop(U%*%c0),0))
  mu=secular(d,b,s)
  return(c(drop(U%*%(b/(d+mu))),mu))
}

# Scalar version for m = 1
row_scalar=function(h,g,s,tau,xold)
{
  if (tau>0)
  {
    h=h+tau
    g=g+tau*xold
  }
  if (g==0) return(c(0,0))
  if (h>0) x0=g/h else x0=Inf
  if (x0^2<=s) return(c(x0,0))
  lo=0
  hi=max(1,abs(g)^(2/3))
  if (h>0) mu=0 else mu=1e-12*hi
  for (it in 1:200)
  {
    f=g^2/(h+mu)^2-s-mu
    if (abs(f)<=1e-15*max(1,s)) break
    if (f>0) lo=mu else hi=mu
    step=mu-f/(-2*g^2/(h+mu)^3-1)
    if (is.finite(step) && step>lo && step<hi) mu=step else mu=0.5*(lo+hi)
    if (hi-lo<=1e-15*max(1,hi)) break
  }
  return(c(g/(h+mu),mu))
}

# PBCD: sweeps over the rows until the gradient norm is <= 1e-9 (or 5 sweeps without decrease)
pbcd=function(S,Lamda,tau=0)
{
  t0=proc.time()[[3]]
  p=nrow(S)
  m=ncol(Lamda)
  s=diag(S)
  val=pbcd_value(S,Lamda)
  stall=0
  for (sweep in 1:100000)
  {
    G=crossprod(Lamda)
    for (i in 1:p)
    {
      li=Lamda[i,]
      gi=drop(crossprod(Lamda,S[,i]))-s[i]*li
      if (m==1) x=row_scalar(G[1,1]-li^2,gi,s[i],tau,li)[1] else x=row_update(G-tcrossprod(li),gi,s[i],tau,li)[1:m]
      G=G-tcrossprod(li)+tcrossprod(x)
      Lamda[i,]=x
    }
    new=pbcd_value(S,Lamda)
    gn=sqrt(sum(pbcd_grad(S,Lamda)^2))
    if (new<val) stall=0 else stall=stall+1
    val=new
    if (gn<=1e-9 || stall>=5) break
  }
  return(list(C=tcrossprod(Lamda),value=val,iter=sweep,conv=gn<=1e-9,sec=proc.time()[[3]]-t0))
}

# L-BFGS-B over psi with one eigendecomposition per point, full or partial (Lanczos)
E3conc_eig=function(S,psi0,m,partial)
{
  p=nrow(S)
  last=NULL
  out=NULL
  ev_all=function(psi)
  {
    if (!is.null(last) && identical(psi,last)) return(out)
    A=S-diag(psi,p)
    if (partial && m<p-1)
    {
      ev=eigs_sym(A,k=m,which="LA",opts=list(tol=1e-13,maxitr=5000,ncv=min(p,max(2*m+1,20))))
      if (length(ev$values)<m) ev=eigen(A,symmetric=TRUE)
      o=order(ev$values,decreasing=TRUE)[1:m]
      val=ev$values[o]
      vec=ev$vectors[,o,drop=FALSE]
    } else {
      ev=eigen(A,symmetric=TRUE)
      val=ev$values[1:m]
      vec=ev$vectors[,1:m,drop=FALSE]
    }
    Lamda=vec%*%diag(sqrt(pmax(val,0)),nrow=m)
    R=A-tcrossprod(Lamda)
    last<<-psi
    out<<-list(value=0.5*sum(R^2),grad=-diag(R),Lamda=Lamda)
    return(out)
  }
  opt=optim(psi0,function(x) ev_all(x)$value,function(x) ev_all(x)$grad,method="L-BFGS-B",
            lower=rep(0,p),upper=diag(S),control=list(maxit=maxit_lbfgs,factr=factr,pgtol=0))
  C=tcrossprod(ev_all(opt$par)$Lamda)
  return(list(C=C,value=Qval(S,C)+Pen(S,C),conv=opt$convergence==0,iter=opt$counts[["function"]]))
}

# Uniquenesses implied by a loading start, inside (0, S_ii)
psi_L=function(S,Lamda)
{
  s=diag(S)
  return(pmin(pmax(s-rowSums(Lamda^2),1e-3*s),0.999*s))
}

meths=c("pbcd","pbcd_prox","e3_conc","e3_conc_one","e3_conc_lanczos","e3_direct","psych_uls")

# Designs: p <= 120 (50 matrices per cell) and p >= 240 (10 matrices per cell, n = 2p)
d1=expand.grid(type=c("T0","T1","D","T2"),n=c(100,500),m=c(1,3,5),p=c(15,30,60,120),stringsAsFactors=FALSE)
d1$R=50
d2=expand.grid(type=c("T0","D","T2"),m=c(1,3,5),p=c(240,480,960),stringsAsFactors=FALSE)
d2$n=2*d2$p
d2$R=10
design=rbind(d1[,c("type","n","m","p","R")],d2[,c("type","n","m","p","R")])
ncell=nrow(design)

# One row per cell, matrix, method and start: cell, rep, method, start, P, gradient norm,
# seconds, iterations, converged, r_max
Res=matrix(NA,sum(design$R)*nst*length(meths),10)
row=0
record=function(cell,rep,meth,j,S,m,C,sec,iter,conv)
{
  Lamda=lamda_C(C,m)
  row<<-row+1
  Res[row,]<<-c(cell,rep,match(meth,meths),j,pbcd_value(S,Lamda),sqrt(sum(pbcd_grad(S,Lamda)^2)),sec,iter,conv,
                max(pmax(diag(C)-diag(S),0)))
}

for (e in 1:ncell)
{
  m=design$m[e]
  p=design$p[e]
  n=design$n[e]
  large=p>=240
  if (design$type[e]=="D") pop=build_pop("A",m,p/m,"high","M2s",NA) else pop=build_pop("A",m,p/m,"moderate",design$type[e],NA)
  R=chol(pop$Sigma)
  for (rep in 1:design$R[e])
  {
    S=draw_S(R,n)
    for (j in 1:nst)
    {
      Lamda0=start_L(S,m,j)
      f=pbcd(S,Lamda0)
      record(e,rep,"pbcd",j,S,m,f$C,f$sec,f$iter,f$conv)
      if (!large)
      {
        f=pbcd(S,Lamda0,tau=1e-3)
        record(e,rep,"pbcd_prox",j,S,m,f$C,f$sec,f$iter,f$conv)
      }
      sec=system.time(f<-E3conc_once(S,psi_L(S,Lamda0),m))[[3]]
      record(e,rep,"e3_conc",j,S,m,f$C,sec,f$iter,f$conv)
      if (large)
      {
        sec=system.time(f<-E3conc_eig(S,psi_L(S,Lamda0),m,FALSE))[[3]]
        record(e,rep,"e3_conc_one",j,S,m,f$C,sec,f$iter,f$conv)
        sec=system.time(f<-E3conc_eig(S,psi_L(S,Lamda0),m,TRUE))[[3]]
        record(e,rep,"e3_conc_lanczos",j,S,m,f$C,sec,f$iter,f$conv)
      } else {
        sec=system.time(f<-E3dir_once(S,Lamda0))[[3]]
        record(e,rep,"e3_direct",j,S,m,f$C,sec,f$iter,f$conv)
      }
    }
    sec=system.time(fa<-suppressMessages(suppressWarnings(fa(r=S,nfactors=m,fm="uls",rotate="none",warnings=FALSE))))[[3]]
    record(e,rep,"psych_uls",1,S,m,tcrossprod(unclass(fa$loadings)[,1:m,drop=FALSE]),sec,NA,1)
  }
}
Res=Res[1:row,]
colnames(Res)=c("cell","rep","method","start","P","gnorm","sec","iter","conv","rmax")

# Per matrix and method: best value over the starts, total time for the five starts
key=paste(Res[,"cell"],Res[,"rep"])
best=tapply(Res[,"P"],key,min)
gap=(Res[,"P"]-best[key])/pmax(1,best[key])
km=paste(key,Res[,"method"])
Per=data.frame(cell=tapply(Res[,"cell"],km,"[",1),method=tapply(Res[,"method"],km,"[",1),
               sec=tapply(Res[,"sec"],km,sum),reached=tapply(gap<=1e-8,km,any),
               gnorm=tapply(Res[,"gnorm"],km,function(x) x[1]),iter=tapply(Res[,"iter"],km,median),
               boundary=tapply(Res[,"rmax"],km,max)>1e-8)
Per$p=design$p[Per$cell]
Per$m=design$m[Per$cell]

# Median time and share of matrices reaching the best value, by method, p and m
Tm=tapply(Per$sec,list(Per$method,Per$p,Per$m),median)
Sh=tapply(Per$reached,list(Per$method,Per$p,Per$m),mean)
dimnames(Tm)[[1]]=meths[as.numeric(dimnames(Tm)[[1]])]
dimnames(Sh)[[1]]=dimnames(Tm)[[1]]
round(Tm,3)
round(Sh,3)
# Matrices in which each method reached the best value
tapply(Per$reached,list(meths[Per$method],Per$p>=240),mean)

# Table S1: p >= 240
TabS=NULL
for (pp in c(240,480,960))
{
  for (mm in c(1,3,5))
  {
    x=Tm[c("pbcd","e3_conc","e3_conc_one","e3_conc_lanczos"),as.character(pp),as.character(mm)]
    TabS=rbind(TabS,c(p=pp,m=mm,x,ratio=x[1]/x[4]))
  }
}
round(TabS,3)

# Figure S1
cols=c(pbcd="#D55E00",pbcd_prox="#E69F00",e3_conc="#0072B2",e3_conc_one="#56B4E9",e3_conc_lanczos="#CC79A7",
       e3_direct="#009E73",psych_uls="grey50")
pchs=c(19,19,15,15,25,17,18)
ltys=c(1,2,1,2,1,1,3)
labs=c("PBCD","PBCD, proximal (tau = 1e-3)","L-BFGS-B over psi, two full eig.","L-BFGS-B, one full eig.",
       "L-BFGS-B, Lanczos","BFGS over Lambda","psych::fa (uls), 1 start")
pv=c(15,30,60,120,240,480,960)
pdf("FigS1.pdf",width=7.5,height=6.5)
layout(matrix(c(1,1,1,2:7),3,byrow=TRUE),heights=c(0.35,1,1))
par(mar=c(0,0,0,0))
plot.new()
legend("center",labs,col=cols,pch=pchs,lty=ltys,pt.bg=cols,ncol=3,bty="n",cex=0.8)
par(mar=c(4,4.5,2,1),mgp=c(2.6,0.7,0))
for (row2 in 1:2)
{
  for (mm in c(1,3,5))
  {
    if (row2==1)
    {
      plot(NA,xlim=c(15,960),ylim=c(0.005,150),log="xy",xaxt="n",yaxt="n",xlab="",
           ylab=ifelse(mm==1,"Median time, 5 starts (s)",""),main=paste("m =",mm),font.main=1,cex.main=1)
      axis(2,at=10^(-2:2),labels=parse(text=paste0("10^",-2:2)),las=1)
    } else {
      plot(NA,xlim=c(15,960),ylim=c(0,1),log="x",xaxt="n",xlab="Number of variables p",
           ylab=ifelse(mm==1,"Share reaching the best value",""),las=1)
    }
    axis(1,at=c(15,60,240,960))
    for (k in 1:7)
    {
      if (row2==1) y=Tm[meths[k],,as.character(mm)] else y=Sh[meths[k],,as.character(mm)]
      ok=!is.na(y)
      lines(pv[ok],y[ok],type="b",col=cols[k],pch=pchs[k],bg=cols[k],lty=ltys[k],lwd=1.5)
    }
  }
}
dev.off()
