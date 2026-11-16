# Bootstrap.R
# Parametric bootstrap test of the residual (Sect. 6.6, Fig. 4). For a sample with an E3 Heywood
# case, resamples of the same n are drawn from the fitted admissible model cov2cor(C + diag(psi)),
# which lies on the boundary; p = (1 + #{r* >= r_obs}) / (B + 1). Samples without a Heywood case
# are not tested and count as non-rejections. Each cell (population and sample size) is saved in
# Bootstrap_cells/ and has its own seed, so the run can be split as in Simulation.R.

source("Estimators.R")

set.seed(123)
B_boot=199
nval=c(100,500)

# Populations: on the boundary (null), interior sampling, structural
pops=data.frame(block=c(rep("A",9),"B"),
                m=c(1,3,1,3,1,1,3,1,1,1),
                k=c(5,3,3,3,5,5,3,5,3,6),
                pattern=c("moderate","moderate","high","high","high","high","moderate","high","high","overshoot"),
                truth=c("T2","T2","T0","T0","M2","T1","T1","M2s","M2s","B"),
                lambda1=c(rep(NA,9),1.2),
                role=c("null","null",rep("interior",5),rep("structural",3)),
                nbsim=c(rep(1000,7),rep(300,3)),stringsAsFactors=FALSE)
npop=nrow(pops)
labels=c("T2, m=1, k=5","T2, m=3, k=3","T0, m=1, k=3, high","T0, m=3, k=3, high","M2, m=1, k=5, high",
         "T1, m=1, k=5, high","T1, m=3, k=3","M2s*, m=1, k=5, high","M2s*, m=1, k=3, high","B, k=6, lambda1 = 1.2")

# Population residual of E3 (positive only for the structural populations)
r_pop=vector()
for (q in 1:npop)
{
  pop=build_pop(pops$block[q],pops$m[q],pops$k[q],pops$pattern[q],pops$truth[q],pops$lambda1[q])
  r_pop[q]=fit_E3(pop$Sigma,pops$m[q],5)$rmax
}
cbind(pops[,c("truth","m","k","role")],r_pop)

# A with A A' = M, for a positive semidefinite M
psd_sqrt=function(M)
{
  ev=eigen(M,symmetric=TRUE)
  return(ev$vectors%*%diag(sqrt(pmax(ev$values,0)),nrow(M)))
}

# Cells: population q and sample size nval[e], numbered cl = 2 (q - 1) + e, with seed 123 + cl.
# "Rscript Bootstrap.R 1 10" runs cells 1 to 10 and "Rscript Bootstrap.R 2 20 6" every sixth cell
# from the second; finished cells are skipped.
args=commandArgs(trailingOnly=TRUE)
ncell=2*npop
step=1
if (length(args)==3) step=as.numeric(args[3])
if (length(args)>=2) cells=seq(as.numeric(args[1]),as.numeric(args[2]),by=step) else cells=1:ncell
dir.create("Bootstrap_cells",showWarnings=FALSE)
cell_file=sprintf("Bootstrap_cells/cell_%02d.rds",1:ncell)
for (cl in cells)
{
  if (file.exists(cell_file[cl])) next
  set.seed(123+cl)
  q=(cl+1)%/%2
  e=cl-2*(q-1)
  pop=build_pop(pops$block[q],pops$m[q],pops$k[q],pops$pattern[q],pops$truth[q],pops$lambda1[q])
  R=chol(pop$Sigma)
  m=pops$m[q]
  p=nrow(pop$Sigma)
  n=nval[e]
  heyw=rep(NA,pops$nbsim[q])
  pval=rep(NA,pops$nbsim[q])
  for (sim in 1:pops$nbsim[q])
  {
    S=draw_S(R,n)
    f=fit_E3(S,m,3,"conc")
    heyw[sim]=f$heywood
    if (isTRUE(f$heywood))
    {
      A=psd_sqrt(cov2cor(f$C+diag(f$psi,p)))
      rstar=vector()
      for (b in 1:B_boot)
      {
        X=matrix(rnorm(n*p),n)%*%t(A)
        rstar[b]=fit_E3(cor(X),m,3,"conc")$rmax
      }
      ok=!is.na(rstar)
      pval[sim]=(1+sum(rstar[ok]>=f$rmax))/(sum(ok)+1)
    }
  }
  saveRDS(list(heyw=heyw,pval=pval),cell_file[cl])
}

# Combine the cells
stopifnot("run the remaining cells first"=all(file.exists(cell_file)))
Heyw=array(NA,c(npop,2,1000))
Pval=array(NA,c(npop,2,1000))
for (cl in 1:ncell)
{
  q=(cl+1)%/%2
  e=cl-2*(q-1)
  x=readRDS(cell_file[cl])
  Heyw[q,e,1:length(x$heyw)]=x$heyw
  Pval[q,e,1:length(x$pval)]=x$pval
}

# Rejection at the 5% level over all samples, with exact (Clopper-Pearson) 95% intervals
Rej=array(NA,c(npop,2,3))
Hrate=matrix(NA,npop,2)
for (q in 1:npop)
{
  for (e in 1:2)
  {
    N=pops$nbsim[q]
    x=sum(Pval[q,e,1:N]<=0.05,na.rm=TRUE)
    Rej[q,e,]=c(x/N,binom.test(x,N)$conf.int)
    Hrate[q,e]=mean(Heyw[q,e,1:N])
  }
}
data.frame(pops[,c("truth","m","k","role")],heywood_100=Hrate[,1],heywood_500=Hrate[,2],
           reject_100=Rej[,1,1],reject_500=Rej[,2,1])

# Pooled rates on the boundary (size), overall and by number of factors
pooled=function(qs)
{
  x=0
  N=0
  for (q in qs)
  {
    x=x+sum(Pval[q,,1:pops$nbsim[q]]<=0.05,na.rm=TRUE)
    N=N+2*pops$nbsim[q]
  }
  return(c(rate=x/N,binom.test(x,N)$conf.int))
}
pooled(1:2)
pooled(1)
pooled(2)
# Rejection among samples with a Heywood case on the boundary
sum(Pval[1:2,,]<=0.05,na.rm=TRUE)/sum(Heyw[1:2,,],na.rm=TRUE)
# Interior sampling populations by n
c(sum(Pval[3:7,1,]<=0.05,na.rm=TRUE),sum(Pval[3:7,2,]<=0.05,na.rm=TRUE))/5000

# Fig. 4: horizontal axis linear below 0.1 and logarithmic above
tf=function(x)
{
  return(ifelse(x<=0.1,x/0.1,1+log10(pmax(x,0.1)/0.1)))
}
ypos=c(13,12,10,9,8,7,6,4,3,2)
ord=c(1,2,5,3,4,6,7,10,9,8)
cols=c(null="#D55E00",interior="#0072B2",structural="black")
cairo_pdf("Fig4.pdf",width=7.5,height=4)
par(mar=c(7,10,1,1),mgp=c(2.4,0.7,0))
plot(NA,xlim=c(0,2),ylim=c(1.2,13.6),axes=FALSE,xlab="Rejection rate at the 5% level over all samples (exact 95% CI)",ylab="")
abline(v=tf(c(0.1,0.2,0.5,1)),col="grey90")
abline(v=tf(0.05),lty=2,col="grey40")
axis(1,at=tf(c(0,0.05,0.1,0.2,0.5,1)),labels=c("0","0.05","0.1","0.2","0.5","1"))
labs=as.expression(labels[ord])
labs[[8]]=quote("B, k=6, "*lambda[1]==1.2)
axis(2,at=ypos,labels=labs,las=1,tick=FALSE,cex.axis=0.8)
for (j in 1:npop)
{
  q=ord[j]
  cl=cols[pops$role[q]]
  for (e in 1:2)
  {
    y=ypos[j]+c(0.15,-0.15)[e]
    arrows(tf(Rej[q,e,2]),y,tf(Rej[q,e,3]),y,angle=90,code=3,length=0.03,col=cl)
    points(tf(Rej[q,e,1]),y,pch=c(19,22)[e],col=cl,bg="white")
  }
}
par(xpd=NA)
legend("bottom",inset=c(0,-0.48),c("On the boundary (H0 true)","Interior sampling populations","Structural populations","n = 100","n = 500"),
       col=c(cols,"grey30","grey30"),pch=c(19,19,19,19,22),ncol=3,bty="n",cex=0.75)
dev.off()
