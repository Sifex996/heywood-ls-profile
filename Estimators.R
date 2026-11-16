# Estimators.R
# Populations and the four criteria compared in the paper. Sourced by the other scripts.
#   E1 unrestricted MINRES, E2 MINRES with C_ii <= S_ii, E3 ULS profile (psi >= 0),
#   E4 Gaussian ML with uniquenesses bounded below by 0.005 (factanal)

# Controls
maxit_bfgs=5000
reltol=1e-16
maxit_bcd=20000
tol_bcd=1e-11
maxit_lbfgs=5000
factr=10
agree_tol=1e-8
tol_bnd=1e-6
ml_lower=0.005

# Populations
build_pop=function(block,m,k,pattern,truth,lambda1)
{
  if (block=="B")
  {
    if (k==3) b=0.25 else b=0.20
    a=lambda1*sqrt(b)
    Sigma=matrix(b,k,k)
    Sigma[1,]=a
    Sigma[,1]=a
    diag(Sigma)=1
    return(list(Sigma=Sigma,psi=NULL))
  }
  if (truth=="M1") mp=m+1 else mp=m
  p=mp*k
  Lamda=matrix(0,p,mp)
  for (f in 1:mp)
  {
    Lamda[((f-1)*k+1):(f*k),f]=seq(0.5,0.7,length.out=k)
  }
  if (pattern=="high") Lamda[1,1]=0.94
  if (truth=="T1") Lamda[1,1]=sqrt(0.98)
  if (truth=="T2") Lamda[1,1]=1
  Phi=matrix(0.3,mp,mp)
  diag(Phi)=1
  K=Lamda%*%Phi%*%t(Lamda)
  psi=1-diag(K)
  Theta=diag(psi,p)
  rho=0
  if (truth=="M2") rho=0.4
  if (truth=="M2s") rho=0.6
  Theta[1,2]=rho*sqrt(psi[1]*psi[2])
  Theta[2,1]=Theta[1,2]
  return(list(Sigma=K+Theta,psi=psi))
}

# Sample correlation matrix of n normal observations (R = chol(Sigma))
draw_S=function(R,n)
{
  Z=matrix(rnorm(n*nrow(R)),n,nrow(R))
  return(cor(Z%*%R))
}

# Criteria
off=function(M)
{
  diag(M)=0
  return(M)
}

Qval=function(S,C)
{
  return(0.5*sum(off(S-C)^2))
}

Pen=function(S,C)
{
  return(0.5*sum(pmax(diag(C)-diag(S),0)^2))
}

# Multiplier of the cap at a capped fit C (zero on the rows where the cap does not bind)
mu_cap=function(S,C)
{
  s=diag(S)
  return(rowSums(off(S-C)*off(C))/s*(diag(C)>=s*(1-tol_bnd)))
}

# Normal-theory covariance matrix of the distinct elements of S (lower triangle, by columns)
gamma_nt=function(S)
{
  idx=which(lower.tri(S,diag=TRUE),arr.ind=TRUE)
  i=idx[,1]
  j=idx[,2]
  return(S[i,i]*S[j,j]+S[i,j]*S[j,i])
}

# Wald statistic of the smallest uniqueness of unrestricted ULS: estimate over standard error,
# from lavaan with normal-theory standard errors. NA when the fit does not converge.
wald_uls=function(S,m,n)
{
  p=nrow(S)
  dimnames(S)=list(paste0("x",1:p),paste0("x",1:p))
  fit=tryCatch(suppressWarnings(lavaan::efa(sample.cov=S,sample.nobs=n,nfactors=m,rotation="none",
                                            estimator="ULS",NACOV=gamma_nt(S),bounds="none",output="lavaan")),
               error=function(e) NULL)
  if (is.null(fit) || !isTRUE(lavaan::lavInspect(fit,"converged"))) return(NA)
  pe=tryCatch(lavaan::parameterEstimates(fit),error=function(e) NULL)
  if (is.null(pe)) return(NA)
  u=pe[pe$op=="~~" & pe$lhs==pe$rhs & substr(pe$lhs,1,1)=="x",]
  i=which.min(u$est)
  z=u$est[i]/u$se[i]
  if (!is.finite(z)) return(NA)
  return(z)
}

# Best rank-m approximation of S - diag(psi)
lamda_psi=function(S,psi,m)
{
  ev=eigen(S-diag(psi,nrow(S)),symmetric=TRUE)
  Lamda=ev$vectors[,1:m,drop=FALSE]%*%diag(sqrt(pmax(ev$values[1:m],0)),nrow=m)
  return(Lamda)
}

# Loadings of a fitted C (orthogonal factors)
lamda_C=function(C,m)
{
  ev=eigen(C,symmetric=TRUE)
  return(ev$vectors[,1:m,drop=FALSE]%*%diag(sqrt(pmax(ev$values[1:m],0)),nrow=m))
}

# Starting values: SMC principal axes (j = 1), then random perturbations
smc=function(S)
{
  s=diag(S)
  inv=tryCatch(solve(S),error=function(e) NULL)
  if (is.null(inv)) psi=0.5*s else psi=1/diag(inv)
  return(pmin(pmax(psi,1e-3*s),0.999*s))
}

start_L=function(S,m,j)
{
  Lamda=lamda_psi(S,smc(S),m)
  if (j>1)
  {
    noise=matrix(rnorm(length(Lamda),sd=0.25),nrow(Lamda),m)
    Lamda=Lamda+sqrt(diag(S))*noise
  }
  return(Lamda)
}

start_psi=function(S,j)
{
  if (j==1) return(smc(S))
  return(runif(nrow(S),0,diag(S)))
}

# min ||s - A x||^2 subject to ||x||^2 <= delta, with H = A'A and g = A's
ls_ball=function(H,g,delta)
{
  m=length(g)
  if (all(g==0)) return(numeric(m))
  eg=eigen(H,symmetric=TRUE)
  U=eg$vectors
  d=pmax(eg$values,0)
  b=drop(crossprod(U,g))
  pos=d>max(d)*1e-12
  c0=numeric(m)
  c0[pos]=b[pos]/d[pos]
  if (!is.finite(delta) || sum(c0^2)<=delta) return(drop(U%*%c0))
  sec=function(mu)
  {
    return(sum(b^2/(d+mu)^2)-delta)
  }
  hi=sqrt(sum(b^2)/delta)
  lo=hi*1e-12
  if (sec(lo)<=0) mu=lo else mu=uniroot(sec,c(lo,hi),tol=1e-15*hi,maxiter=1000)$root
  x=drop(U%*%(b/(d+mu)))
  return(x*sqrt(delta/sum(x^2)))
}

# Gauss-Seidel over loading rows (Harman-Jones); caps = Inf: E1, caps = diag(S): E2
minres_bcd=function(S,Lamda,caps)
{
  p=nrow(S)
  nr=rowSums(Lamda^2)
  Lamda=Lamda*ifelse(is.finite(caps) & nr>caps,sqrt(caps/nr),1)
  conv=FALSE
  it=0
  for (it in 1:maxit_bcd)
  {
    Lold=Lamda
    for (i in 1:p)
    {
      A=Lamda[-i,,drop=FALSE]
      Lamda[i,]=ls_ball(crossprod(A),drop(crossprod(A,S[i,-i])),caps[i])
    }
    if (max(abs(Lamda-Lold))<tol_bcd)
    {
      conv=TRUE
      break
    }
  }
  return(list(Lamda=Lamda,iter=it,conv=conv))
}

# BFGS over the loadings on P_w(L) = Q(L) + w/2 sum_i (C_ii - S_ii)_+^2
fit_Pw=function(S,Lamda0,w,rtol=reltol)
{
  p=nrow(S)
  m=ncol(Lamda0)
  s=diag(S)
  crit=function(x)
  {
    Lamda=matrix(x,p,m)
    C=tcrossprod(Lamda)
    R=off(S-C)
    val=0.5*sum(R^2)
    if (w>0)
    {
      e=pmax(diag(C)-s,0)
      val=val+0.5*w*sum(e^2)
    }
    return(val)
  }
  grad=function(x)
  {
    Lamda=matrix(x,p,m)
    C=tcrossprod(Lamda)
    R=off(S-C)
    G=-2*R%*%Lamda
    if (w>0)
    {
      e=pmax(diag(C)-s,0)
      G=G+2*w*e*Lamda
    }
    return(as.vector(G))
  }
  opt=optim(as.vector(Lamda0),crit,grad,method="BFGS",control=list(maxit=maxit_bfgs,reltol=rtol))
  Lamda=matrix(opt$par,p,m)
  return(list(Lamda=Lamda,value=opt$value,conv=opt$convergence==0,iter=opt$counts[["function"]],
              gnorm=sqrt(sum(grad(opt$par)^2))))
}

# One start of each criterion
E1_once=function(S,Lamda0)
{
  b=fit_Pw(S,Lamda0,0)
  pol=minres_bcd(S,b$Lamda,rep(Inf,nrow(S)))
  C=tcrossprod(pol$Lamda)
  return(list(C=C,value=Qval(S,C),conv=b$conv || pol$conv,iter=b$iter+pol$iter))
}

E2_once=function(S,Lamda0)
{
  r=minres_bcd(S,Lamda0,diag(S))
  C=tcrossprod(r$Lamda)
  return(list(C=C,value=Qval(S,C),conv=r$conv,iter=r$iter))
}

# E3 over psi in [0, S_ii], loadings concentrated out (L-BFGS-B)
E3conc_once=function(S,psi0,m)
{
  s=diag(S)
  crit=function(psi)
  {
    Lamda=lamda_psi(S,psi,m)
    R=S-diag(psi,nrow(S))-tcrossprod(Lamda)
    return(0.5*sum(R^2))
  }
  grad=function(psi)
  {
    Lamda=lamda_psi(S,psi,m)
    R=S-diag(psi,nrow(S))-tcrossprod(Lamda)
    return(-diag(R))
  }
  opt=optim(psi0,crit,grad,method="L-BFGS-B",lower=rep(0,length(s)),upper=s,
            control=list(maxit=maxit_lbfgs,factr=factr,pgtol=0))
  C=tcrossprod(lamda_psi(S,opt$par,m))
  return(list(C=C,value=Qval(S,C)+Pen(S,C),conv=opt$convergence==0,
              iter=opt$counts[["function"]],psi_opt=opt$par))
}

# E3 over the loadings on the profile (BFGS)
E3dir_once=function(S,Lamda0)
{
  b=fit_Pw(S,Lamda0,1)
  C=tcrossprod(b$Lamda)
  return(list(C=C,value=Qval(S,C)+Pen(S,C),conv=b$conv,iter=b$iter))
}

# Keep the best start; n_agree counts the starts within agree_tol of it
keep_best=function(fits)
{
  if (length(fits)==0) return(NULL)
  fits=fits[!sapply(fits,is.null)]
  if (length(fits)==0) return(NULL)
  vals=sapply(fits,function(f) f$value)
  best=fits[[which.min(vals)]]
  best$n_agree=sum(vals<=min(vals)+agree_tol*max(1,abs(min(vals))))
  best$n_ok=length(fits)
  return(best)
}

# Common report: psi as each criterion reports it, Heywood flags, residual r
report=function(S,raw,psi,bnd=NULL)
{
  s=diag(S)
  if (is.null(bnd)) bnd=psi<=tol_bnd*s
  r=pmax(diag(raw$C)-s,0)
  return(list(C=raw$C,psi=psi,value=raw$value,Q=Qval(S,raw$C),pen=Pen(S,raw$C),
              bnd=bnd,heywood=any(bnd),improper=any(psi< -tol_bnd*s),r=r,rmax=max(r),
              negpsi=max(pmax(-psi,0)),conv=isTRUE(raw$conv),iter=raw$iter,
              n_agree=raw$n_agree,n_ok=raw$n_ok,psi_opt=raw$psi_opt))
}

failed=function(S)
{
  p=nrow(S)
  return(list(C=matrix(NA,p,p),psi=rep(NA,p),value=NA,Q=NA,pen=NA,bnd=rep(NA,p),heywood=NA,
              improper=NA,r=rep(NA,p),rmax=NA,negpsi=NA,conv=FALSE,iter=NA,n_agree=NA,n_ok=0))
}

# Estimators (nst starting values)
fit_E1=function(S,m,nst)
{
  fits=list()
  for (j in 1:nst)
  {
    Lamda0=start_L(S,m,j)
    fits[[j]]=tryCatch(E1_once(S,Lamda0),error=function(e) NULL)
  }
  raw=keep_best(fits)
  if (is.null(raw)) return(failed(S))
  return(report(S,raw,diag(S)-diag(raw$C)))
}

fit_E2=function(S,m,nst)
{
  fits=list()
  for (j in 1:nst)
  {
    Lamda0=start_L(S,m,j)
    fits[[j]]=tryCatch(E2_once(S,Lamda0),error=function(e) NULL)
  }
  raw=keep_best(fits)
  if (is.null(raw)) return(failed(S))
  return(report(S,raw,pmax(diag(S)-diag(raw$C),0)))
}

# variant "best": smaller value of the two algorithms; "conc" or "direct": one only
fit_E3=function(S,m,nst,variant="best")
{
  conc=NULL
  dir=NULL
  if (variant!="direct")
  {
    fits=list()
    for (j in 1:nst)
    {
      psi0=start_psi(S,j)
      fits[[j]]=tryCatch(E3conc_once(S,psi0,m),error=function(e) NULL)
    }
    conc=keep_best(fits)
  }
  if (variant!="conc")
  {
    fits=list()
    for (j in 1:nst)
    {
      Lamda0=start_L(S,m,j)
      fits[[j]]=tryCatch(E3dir_once(S,Lamda0),error=function(e) NULL)
    }
    dir=keep_best(fits)
  }
  if (is.null(conc) && is.null(dir)) return(failed(S))
  raw=conc
  if (is.null(conc)) raw=dir
  if (!is.null(conc) && !is.null(dir) && dir$value<conc$value) raw=dir
  return(report(S,raw,pmax(diag(S)-diag(raw$C),0)))
}

fit_E4=function(S,m,n)
{
  s=diag(S)
  fa=tryCatch(suppressWarnings(factanal(factors=m,covmat=S,n.obs=n,rotation="none",lower=ml_lower)),
              error=function(e) NULL)
  if (is.null(fa)) return(failed(S))
  Lamda=unclass(fa$loadings)[,1:m,drop=FALSE]*sqrt(s)
  u=fa$uniquenesses
  raw=list(C=tcrossprod(Lamda),value=unname(fa$criteria[["objective"]]),conv=isTRUE(fa$converged),
           iter=NA,n_agree=NA,n_ok=NA)
  return(report(S,raw,u*s,bnd=u<=ml_lower+1e-4))
}
