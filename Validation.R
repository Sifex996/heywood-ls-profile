# Validation.R
# Checks of the estimators before the simulation (Sect. 6.1-6.2): populations, the 3 x 3
# boundary example, agreement of the two E3 algorithms, E1 against psych::fa, interior E3 = E1,
# the profile identity, and one fit of every population. Every gate should show pass = cases.
# The last block checks whether the E3 fit of each population matrix is unique.

library("psych")
source("Estimators.R")

set.seed(123)
gate=data.frame(gate=character(),cases=numeric(),pass=numeric(),max_disc=numeric())

# G0: populations positive definite with psi > 0; E1 and E3 fit T0, T1 exactly; Block B structural
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
npop=nrow(pops)
mineig=vector()
minpsi=vector()
r_pop=vector()
val_exact=vector()
bnd_exact=vector()
for (q in 1:npop)
{
  pop=build_pop(pops$block[q],pops$m[q],pops$k[q],pops$pattern[q],pops$truth[q],pops$lambda1[q])
  mineig[q]=min(eigen(pop$Sigma,symmetric=TRUE,only.values=TRUE)$values)
  if (is.null(pop$psi)) minpsi[q]=NA else minpsi[q]=min(pop$psi)
  f1=fit_E1(pop$Sigma,pops$m[q],5)
  f3=fit_E3(pop$Sigma,pops$m[q],5)
  r_pop[q]=f3$rmax
  val_exact[q]=max(f1$value,f3$value)
  bnd_exact[q]=f1$heywood || f3$heywood
}
ok=mineig>0 & (is.na(minpsi) | minpsi>0)
gate[1,]=list("G0a",npop,sum(ok),-min(c(mineig,0)))
corr=pops$truth %in% c("T0","T1")
gate[2,]=list("G0b",sum(corr),sum(val_exact[corr]<1e-9 & !bnd_exact[corr]),max(val_exact[corr]))
gate[3,]=list("G0c",sum(pops$block=="B"),sum(r_pop[pops$block=="B"]>1e-6),0)
sum(r_pop>1e-6)

# G1: 3 x 3 example of Sect. 6.1 (closed-form reference values)
a=0.60
b=0.25
S3=build_pop("B",1,3,"overshoot","B",1.2)$Sigma
t_cap=uniroot(function(t) t^3+(1-b)*t-a,c(sqrt(b),a),tol=1e-14)$root
x2=function(q)
{
  return((1+sqrt(1+8*q*(q-b)))/2)
}
q_soft=uniroot(function(q) sqrt(q)*(sqrt(x2(q))+(q-b)/sqrt(x2(q)))-a,c(b,a^2),tol=1e-14)$root
C1=tcrossprod(c(a/sqrt(b),sqrt(b),sqrt(b)))
C2=tcrossprod(c(1,t_cap,t_cap))
C3=tcrossprod(c(sqrt(x2(q_soft)),sqrt(q_soft),sqrt(q_soft)))
M_H=Qval(S3,C2)
M_P=Qval(S3,C3)+Pen(S3,C3)
f1=fit_E1(S3,1,5)
f2=fit_E2(S3,1,5)
f3=fit_E3(S3,1,5)
ok1=c(abs(M_H-0.0072222648)<1e-10 && abs(M_P-0.0064642016)<1e-10,
      abs(f1$value)<1e-12 && max(abs(f1$C-C1))<1e-5,
      abs(f2$value-M_H)<1e-8 && max(abs(f2$C-C2))<1e-5,
      abs(f3$value-M_P)<1e-8 && max(abs(f3$C-C3))<1e-5,
      f2$bnd[1] && f3$bnd[1] && f1$improper,
      abs(f3$rmax-(x2(q_soft)-1))<1e-6)
gate[4,]=list("G1",6,sum(ok1),max(abs(f2$C-C2),abs(f3$C-C3)))
c(f2$value,M_H,f3$value,M_P,f3$rmax)

# Random correlation matrices with communalities in comm
random_S=function(p,m,n,comm)
{
  L=matrix(runif(p*m,-1,1),p,m)
  L=L/sqrt(rowSums(L^2))*sqrt(runif(p,comm[1],comm[2]))
  Sigma=tcrossprod(L)
  diag(Sigma)=1
  return(cor(matrix(rnorm(n*p),n,p)%*%chol(Sigma)))
}

# G2: the two E3 algorithms reach the same value (6 x 6, n = 60, frequent boundary)
S2=list()
d2=vector()
bnd2=vector()
conc=list()
for (i in 1:200)
{
  S2[[i]]=random_S(6,2,60,c(0.3,0.9))
  conc[[i]]=fit_E3(S2[[i]],2,5,"conc")
  dir=fit_E3(S2[[i]],2,5,"direct")
  d2[i]=abs(conc[[i]]$value-dir$value)/max(1,abs(dir$value))
  bnd2[i]=conc[[i]]$heywood
}
gate[5,]=list("G2",200,sum(d2<=1e-7),max(d2))
sum(bnd2)

# G3: E1 against psych::fa(fm = "minres"); G4: interior E3 = E1
S34=list()
q_gap=vector()
f1s=list()
f3s=list()
for (i in 1:50)
{
  S34[[i]]=random_S(8,2,1000,c(0.3,0.6))
  f1s[[i]]=fit_E1(S34[[i]],2,5)
  fa=suppressMessages(suppressWarnings(fa(r=S34[[i]],nfactors=2,n.obs=1000,rotate="none",fm="minres",warnings=FALSE)))
  Lp=unclass(fa$loadings)[,1:2]
  q_gap[i]=f1s[[i]]$Q-Qval(S34[[i]],tcrossprod(Lp))
}
gate[6,]=list("G3",50,sum(q_gap<=1e-8),max(q_gap))
d4=vector()
ok4=vector()
for (i in 1:50)
{
  f3s[[i]]=fit_E3(S34[[i]],2,5)
  if (!f3s[[i]]$heywood)
  {
    dv=abs(f3s[[i]]$value-f1s[[i]]$value)
    dC=max(abs(f3s[[i]]$C-f1s[[i]]$C))
    d4=c(d4,max(dv,dC*1e-3))
    ok4=c(ok4,dv<=1e-7*max(1,f1s[[i]]$value) && dC<=1e-4)
  }
}
gate[7,]=list("G4",length(ok4),sum(ok4),max(d4))

# G5: F_ULS = Q + penalty, and the optimizer's psi equals (S_ii - C_ii)_+
f5=c(conc,f3s)
S5=c(S2,S34)
d5=vector()
for (i in 1:length(f5))
{
  ident=abs(0.5*sum((S5[[i]]-f5[[i]]$C-diag(f5[[i]]$psi))^2)-(f5[[i]]$Q+f5[[i]]$pen))
  if (is.null(f5[[i]]$psi_opt)) kkt=0 else kkt=max(abs(f5[[i]]$psi_opt-f5[[i]]$psi))
  d5[i]=max(ident/1e-12,kkt/1e-5)
}
gate[8,]=list("G5",length(f5),sum(d5<=1),max(d5))

# G6: all four criteria on one sample (n = 200) of every population, with time per sample
Sec=vector()
fail=vector()
for (q in 1:npop)
{
  pop=build_pop(pops$block[q],pops$m[q],pops$k[q],pops$pattern[q],pops$truth[q],pops$lambda1[q])
  S=draw_S(chol(pop$Sigma),200)
  m=pops$m[q]
  if (m==3) nst=10 else nst=5
  t0=proc.time()[[3]]
  v=c(fit_E1(S,m,nst)$value,fit_E2(S,m,nst)$value,fit_E3(S,m,nst)$value,fit_E4(S,m,200)$value)
  Sec[q]=proc.time()[[3]]-t0
  fail[q]=any(is.na(v))
}
gate[9,]=list("G6",npop,sum(!fail),0)

gate
# Projected single-core hours of Simulation.R (5 sample sizes, 1000 samples)
sum(Sec)*5*1000/3600

# Uniqueness of the E3 population fit (Sect. 6.2): 40 starting values for each of the two
# algorithms. at_min counts the starts within 1e-9 of the smallest value, max_diff_C is the
# largest difference between their common covariances and max_r their largest residual.
uniq=NULL
for (q in 1:npop)
{
  pop=build_pop(pops$block[q],pops$m[q],pops$k[q],pops$pattern[q],pops$truth[q],pops$lambda1[q])
  S=pop$Sigma
  m=pops$m[q]
  vals=vector()
  Cs=list()
  for (j in 1:40)
  {
    f=fit_Pw(S,start_L(S,m,j),1,1e-15)
    vals=c(vals,f$value)
    Cs[[length(Cs)+1]]=tcrossprod(f$Lamda)
    f=E3conc_once(S,start_psi(S,j),m)
    vals=c(vals,f$value)
    Cs[[length(Cs)+1]]=f$C
  }
  best=which.min(vals)
  near=which(vals<=min(vals)+1e-9)
  dC=0
  rr=0
  for (j in near)
  {
    dC=max(dC,max(abs(Cs[[j]]-Cs[[best]])))
    rr=max(rr,max(pmax(diag(Cs[[j]])-diag(S),0)))
  }
  uniq=rbind(uniq,data.frame(pops[q,],starts=length(vals),at_min=length(near),max_diff_C=dC,max_r=rr))
}
uniq
# Populations with one common covariance, and how many of them are structural
c(sum(uniq$max_diff_C<=1e-6),sum(uniq$max_diff_C<=1e-6 & r_pop>1e-6))
uniq[uniq$max_diff_C>1e-6,]
