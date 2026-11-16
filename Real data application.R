# Real data application.R
# Sect. 7: Heywood cases in correlation matrices and data sets distributed with R and the packages
# psych, lavaan and MASS. Two screens, Tables 5 and 6 with the parametric bootstrap check (B = 499),
# the spread over starting values, and the fit of psych::fa for comparison.

library("psych")
library("lavaan")
library("MASS")
source("Estimators.R")

set.seed(123)
B_boot=499

# Correlation matrix of the numeric columns (complete cases) of a data set
raw_cor=function(x,cols=NULL)
{
  x=as.data.frame(x)
  x=x[,sapply(x,is.numeric),drop=FALSE]
  if (!is.null(cols)) x=x[,cols,drop=FALSE]
  x=x[complete.cases(x),,drop=FALSE]
  return(list(S=cor(x),n=nrow(x)))
}

# E1-E3 on a correlation matrix, with the reports of Tables 5 and 6
fit3=function(S,m,nst)
{
  vars=colnames(S)
  S=cov2cor((S+t(S))/2)
  dimnames(S)=NULL
  f1=fit_E1(S,m,nst)
  f2=fit_E2(S,m,nst)
  f3=fit_E3(S,m,nst)
  i=which.max(f3$r)
  # multiplier of the cap at the capped fit (active rows only)
  mu=mu_cap(S,f2$C)
  return(list(S=S,f1=f1,f2=f2,f3=f3,row=data.frame(variable=ifelse(is.null(vars),i,vars[i]),
    n_heywood=sum(f3$bnd),heywood_E1=f1$heywood,heywood_E2=f2$heywood,heywood_E3=f3$heywood,
    psi_E1=f1$psi[i],r_max=f3$rmax,mu=max(mu),bound=sqrt(2*max(f2$Q-f3$Q,0)),M0=f1$Q,MP=f3$value,MH=f2$Q,
    comm_diff=max(abs(diag(f3$C)-diag(f2$C))),off_rms=sqrt(mean((f3$C-f2$C)[upper.tri(S)]^2)),
    agree_E3=f3$n_agree)))
}

# Parametric bootstrap check of the residual from the fitted admissible model
boot_p=function(S,f3,n,m)
{
  p=nrow(S)
  ev=eigen(cov2cor(f3$C+diag(f3$psi,p)),symmetric=TRUE)
  A=ev$vectors%*%(sqrt(pmax(ev$values,0))*t(ev$vectors))
  rstar=vector()
  for (b in 1:B_boot)
  {
    X=matrix(rnorm(n*p),n,p)%*%A
    rstar[b]=tryCatch(fit_E3(cor(X),m,3,"conc")$rmax,error=function(e) NA)
  }
  ok=!is.na(rstar)
  return((1+sum(rstar[ok]>=f3$rmax))/(sum(ok)+1))
}

# First screen: ten published correlation matrices ----------------------------------------------
screen1=list(
  list(name="Harman74 (24 tests)",S=Harman74.cor$cov,n=145,m=3:5),
  list(name="ability.cov (6 tests)",S=cov2cor(ability.cov$cov),n=112,m=1:2),
  list(name="Harman.8 (8 physical variables)",S=Harman.8,n=305,m=1:3),
  list(name="Harman.5 (5 socio-economic variables)",S=cor(Harman.5),n=12,m=1:2),
  list(name="Holzinger.9 (9 tests)",S=Holzinger.9,n=145,m=2:4),
  list(name="Thurstone (9 tests)",S=Thurstone,n=213,m=2:4),
  list(name="Thurstone.33 (9 tests)",S=Thurstone.33,n=4175,m=2:4),
  list(name="Bechtoldt.1 (17 tests)",S=Bechtoldt.1,n=212,m=4:6),
  list(name="Reise (16 items)",S=Reise,n=35000,m=3:5),
  list(name="Tucker (9 variables)",S=Tucker,n=710,m=2:3))
Scr1=NULL
for (cs in screen1)
{
  for (m in cs$m)
  {
    p=nrow(cs$S)
    if ((p-m)^2-p-m>=0)
    {
      f=fit3(as.matrix(cs$S),m,10)
      Scr1=rbind(Scr1,data.frame(data=cs$name,p=p,n=cs$n,m=m,f$row))
    }
  }
}
Scr1[,c("data","m","heywood_E1","heywood_E2","heywood_E3","psi_E1","r_max")]
c(fits=nrow(Scr1),heywood=sum(Scr1$heywood_E3))

# Table 5 --------------------------------------------------------------------------------------------
cases5=list(
  list(id="Harman5_m1",S=cor(Harman.5),n=12,m=1),
  list(id="Harman5_m2",S=cor(Harman.5),n=12,m=2),
  list(id="Harman8_m3",S=Harman.8,n=305,m=3),
  list(id="Holzinger9_m4",S=Holzinger.9,n=145,m=4),
  list(id="Thurstone_m4",S=Thurstone,n=213,m=4),
  list(id="Bechtoldt1_m6",S=Bechtoldt.1,n=212,m=6))
Tab5=NULL
for (cs in cases5)
{
  f=fit3(as.matrix(cs$S),cs$m,10)
  Tab5=rbind(Tab5,data.frame(id=cs$id,p=nrow(f$S),n=cs$n,m=cs$m,f$row,boot_p=boot_p(f$S,f$f3,cs$n,cs$m)))
}
Tab5[,c("id","p","n","m","variable","n_heywood","psi_E1","r_max","mu","bound","boot_p")]
Tab5[,c("id","M0","MP","MH","comm_diff","off_rms","agree_E3")]

# Spread over 50 starting values: profile (w = 1) and unrestricted MINRES (w = 0)
Starts=NULL
for (cs in cases5[-1])
{
  S=cov2cor((as.matrix(cs$S)+t(as.matrix(cs$S)))/2)
  dimnames(S)=NULL
  v3=vector()
  r3=vector()
  v1=vector()
  p1=vector()
  for (j in 1:50)
  {
    Lamda0=start_L(S,cs$m,j)
    f=fit_Pw(S,Lamda0,1,1e-15)
    v3[j]=f$value
    r3[j]=max(pmax(rowSums(f$Lamda^2)-diag(S),0))
    f=fit_Pw(S,Lamda0,0,1e-15)
    v1[j]=f$value
    p1[j]=min(diag(S)-rowSums(f$Lamda^2))
  }
  near3=v3<=min(v3)+1e-8
  near1=v1<=min(v1)+1e-8
  Starts=rbind(Starts,data.frame(id=cs$id,E3_at_min=sum(near3),E3_r_lo=min(r3[near3]),E3_r_hi=max(r3[near3]),
                                 E1_at_min=sum(near1),E1_psi_lo=min(p1[near1]),E1_psi_hi=max(p1[near1]),E1_min=min(v1)))
}
Starts

# psych::fa (minres, lower bound 0.005 on the uniquenesses) on the same five fits
for (cs in cases5[-1])
{
  fa=suppressWarnings(suppressMessages(fa(cs$S,nfactors=cs$m,rotate="none",fm="minres")))
  print(c(cs$id,round(min(fa$uniquenesses),5)))
}

# Second screen: 30 further data sets, one to three factors -------------------------------------------
screen2=list(
  "psych Harman.political"=list(S=as.matrix(Harman.political),n=147),
  "psych Dwyer"=list(S=as.matrix(Dwyer),n=1000),
  "psych cattell"=list(S=as.matrix(cattell),n=277),
  "psych Holzinger"=list(S=as.matrix(Holzinger),n=355),
  "psych Bechtoldt.2"=list(S=as.matrix(Bechtoldt.2),n=213),
  "psych Thurstone.9"=list(S=as.matrix(Thurstone.9),n=NA),
  "psych bfi"=raw_cor(bfi,1:25),
  "psych sat.act"=raw_cor(sat.act,c("ACT","SATV","SATQ")),
  "lavaan HolzingerSwineford1939"=raw_cor(HolzingerSwineford1939,paste0("x",1:9)),
  "lavaan PoliticalDemocracy"=raw_cor(PoliticalDemocracy),
  "datasets USJudgeRatings"=raw_cor(USJudgeRatings),
  "datasets mtcars"=raw_cor(mtcars),
  "datasets swiss"=raw_cor(swiss),
  "datasets state.x77"=raw_cor(state.x77),
  "datasets attitude"=raw_cor(attitude),
  "datasets LifeCycleSavings"=raw_cor(LifeCycleSavings),
  "datasets longley"=raw_cor(longley),
  "datasets iris"=raw_cor(iris),
  "datasets trees"=raw_cor(trees),
  "datasets stackloss"=raw_cor(stackloss),
  "datasets rock"=raw_cor(rock),
  "datasets airquality"=raw_cor(airquality,1:4),
  "datasets freeny"=raw_cor(freeny),
  "MASS crabs"=raw_cor(crabs,c("FL","RW","CL","CW","BD")),
  "MASS UScrime"=raw_cor(UScrime),
  "MASS Boston"=raw_cor(Boston),
  "MASS painters"=raw_cor(painters),
  "MASS hills"=raw_cor(hills),
  "MASS cpus"=raw_cor(cpus),
  "MASS fgl"=raw_cor(fgl))
Scr2=NULL
for (nm in names(screen2))
{
  p=nrow(screen2[[nm]]$S)
  for (m in 1:3)
  {
    if ((p-m)^2-p-m>=0)
    {
      f=fit3(screen2[[nm]]$S,m,20)
      Scr2=rbind(Scr2,data.frame(data=nm,p=p,n=screen2[[nm]]$n,m=m,f$row))
    }
  }
}
Scr2[Scr2$heywood_E3,c("data","p","n","m","variable","n_heywood","psi_E1","r_max")]
c(data=length(unique(Scr2$data)),fits=nrow(Scr2),heywood=sum(Scr2$heywood_E3),r_above_0.01=sum(Scr2$r_max>0.01))

# Table 6 --------------------------------------------------------------------------------------------
cases6=list(
  list(id="political_m3",S=as.matrix(Harman.political),n=147,m=3),
  list(id="iris_m1",S=raw_cor(iris)$S,n=150,m=1),
  list(id="fgl_m1",S=raw_cor(fgl)$S,n=214,m=1),
  list(id="rock_m1",S=raw_cor(rock)$S,n=48,m=1),
  list(id="hills_m1",S=raw_cor(hills)$S,n=35,m=1),
  list(id="trees_m1",S=raw_cor(trees)$S,n=31,m=1),
  list(id="longley_m1",S=raw_cor(longley)$S,n=16,m=1))
Tab6=NULL
for (cs in cases6)
{
  f=fit3(cs$S,cs$m,20)
  Tab6=rbind(Tab6,data.frame(id=cs$id,p=nrow(f$S),n=cs$n,m=cs$m,f$row,boot_p=boot_p(f$S,f$f3,cs$n,cs$m)))
}
Tab6[,c("id","p","n","m","variable","n_heywood","psi_E1","r_max","mu","bound","boot_p")]
