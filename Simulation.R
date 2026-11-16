# Simulation.R
# Simulation of Sect. 6.2-6.5: 46 populations, n in {50, 100, 200, 500, 1000},
# nbsim samples per cell, criteria E1-E4 on every sample correlation matrix, and for the
# samples with a Heywood case the multiplier of the cap and a Wald statistic.
# Each design cell is saved in Simulation_cells/ and has its own seed, 123 + cell number, so the
# run can be split: "Rscript Simulation.R 1 40" runs cells 1 to 40, "Rscript Simulation.R 2 230 6"
# runs every sixth cell from the second, and several such sessions can run side by side.
# Finished cells are skipped. Without arguments all 230 cells are run. When
# every cell is there, the results are combined in Simulation_results.rds for
# "Simulation summary.R". The numbers do not depend on how the run is split.

source("Estimators.R")

nbsim=1000
nval=c(50,100,200,500,1000)

# Design: Block A (full factorial) and Block B (overshoot matrices)
A=expand.grid(n=nval,truth=c("T0","T1","M1","M2","M2s"),pattern=c("moderate","high"),k=c(3,5),m=c(1,3),
              stringsAsFactors=FALSE)
A$block="A"
A$lambda1=NA
B=expand.grid(n=nval,lambda1=c(1.1,1.2,1.3),k=c(3,6),stringsAsFactors=FALSE)
B$m=1
B$truth="B"
B$pattern="overshoot"
B$block="B"
keys=c("block","m","k","pattern","truth","lambda1")
design=rbind(A[,c(keys,"n")],B[,c(keys,"n")])
design=design[order(design$block,design$m,design$k,design$pattern,design$truth,design$lambda1,design$n),]
rownames(design)=NULL
ncell=nrow(design)

# Cells of this session
args=commandArgs(trailingOnly=TRUE)
step=1
if (length(args)==3) step=as.numeric(args[3])
if (length(args)>=2) cells=seq(as.numeric(args[1]),as.numeric(args[2]),by=step) else cells=1:ncell

# Populations (cells that differ only in n)
pops=unique(design[,keys])
rownames(pops)=NULL
npop=nrow(pops)
design$pop=match(paste(design$block,design$m,design$k,design$pattern,design$truth,design$lambda1),
                 paste(pops$block,pops$m,pops$k,pops$pattern,pops$truth,pops$lambda1))

# Population fits (n = infinity): label "structural" when the E3 residual is positive
set.seed(123)
r_pop=vector()
negpsi_pop=vector()
mu_pop=vector()
Cpop=list()
for (q in 1:npop)
{
  pop=build_pop(pops$block[q],pops$m[q],pops$k[q],pops$pattern[q],pops$truth[q],pops$lambda1[q])
  f1=fit_E1(pop$Sigma,pops$m[q],5)
  f2=fit_E2(pop$Sigma,pops$m[q],5)
  f3=fit_E3(pop$Sigma,pops$m[q],5)
  r_pop[q]=f3$rmax
  negpsi_pop[q]=f1$negpsi
  mu_pop[q]=max(mu_cap(pop$Sigma,f2$C))
  Cpop[[q]]=list(f1$C,f2$C,f3$C)
}
structural=r_pop>1e-6

# One file per cell, one row per sample
dir.create("Simulation_cells",showWarnings=FALSE)
cell_file=sprintf("Simulation_cells/cell_%03d.rds",1:ncell)
for (e in cells)
{
  if (file.exists(cell_file[e])) next
  set.seed(123+e)
  pop=build_pop(design$block[e],design$m[e],design$k[e],design$pattern[e],design$truth[e],design$lambda1[e])
  R=chol(pop$Sigma)
  m=design$m[e]
  n=design$n[e]
  if (m==3) nst=10 else nst=5
  Hey=matrix(NA,nbsim,4)
  Conv=matrix(NA,nbsim,4)
  Agree=matrix(NA,nbsim,3)
  Err=matrix(NA,nbsim,3)
  Rmse=matrix(NA,nbsim,3)
  Imp=rep(NA,nbsim)
  Neg=rep(NA,nbsim)
  Rmax=rep(NA,nbsim)
  Psimin=rep(NA,nbsim)
  Comm=rep(NA,nbsim)
  Mu=rep(NA,nbsim)
  Wald=rep(NA,nbsim)
  Sec=rep(NA,nbsim)
  for (sim in 1:nbsim)
  {
    S=draw_S(R,n)
    t0=proc.time()[[3]]
    f1=fit_E1(S,m,nst)
    f2=fit_E2(S,m,nst)
    f3=fit_E3(S,m,nst)
    f4=fit_E4(S,m,n)
    Sec[sim]=proc.time()[[3]]-t0
    fits=list(f1,f2,f3,f4)
    for (k in 1:4)
    {
      Hey[sim,k]=fits[[k]]$heywood
      Conv[sim,k]=fits[[k]]$conv
    }
    # Recovery of the population fit of the same criterion
    for (k in 1:3)
    {
      D=fits[[k]]$C-Cpop[[design$pop[e]]][[k]]
      Err[sim,k]=max(abs(diag(D)))
      Rmse[sim,k]=sqrt(mean(D[upper.tri(D)]^2))
      Agree[sim,k]=fits[[k]]$n_agree==fits[[k]]$n_ok
    }
    Imp[sim]=f1$improper
    Neg[sim]=f1$negpsi
    Rmax[sim]=f3$rmax
    Psimin[sim]=min(f3$psi)
    Comm[sim]=max(abs(diag(f3$C)-diag(f2$C)))
    # Multiplier of the cap (zero when no cap binds), and the Wald statistic for Heywood samples
    Mu[sim]=max(mu_cap(S,f2$C))
    if (isTRUE(f3$heywood)) Wald[sim]=wald_uls(S,m,n)
  }
  saveRDS(list(nbsim=nbsim,Hey=Hey,Conv=Conv,Agree=Agree,Err=Err,Rmse=Rmse,Imp=Imp,Neg=Neg,Rmax=Rmax,
               Psimin=Psimin,Comm=Comm,Mu=Mu,Wald=Wald,Sec=Sec),cell_file[e])
}

# Combine the cells, one row per cell and one column per sample
if (all(file.exists(cell_file)))
{
  Hey=array(NA,c(ncell,nbsim,4))
  Conv=array(NA,c(ncell,nbsim,4))
  Agree=array(NA,c(ncell,nbsim,3))
  Err=array(NA,c(ncell,nbsim,3))
  Rmse=array(NA,c(ncell,nbsim,3))
  Imp=matrix(NA,ncell,nbsim)
  Neg=matrix(NA,ncell,nbsim)
  Rmax=matrix(NA,ncell,nbsim)
  Psimin=matrix(NA,ncell,nbsim)
  Comm=matrix(NA,ncell,nbsim)
  Mu=matrix(NA,ncell,nbsim)
  Wald=matrix(NA,ncell,nbsim)
  Sec=matrix(NA,ncell,nbsim)
  for (e in 1:ncell)
  {
    x=readRDS(cell_file[e])
    if (x$nbsim!=nbsim) stop("cell file with another number of samples: ",cell_file[e])
    Hey[e,,]=x$Hey
    Conv[e,,]=x$Conv
    Agree[e,,]=x$Agree
    Err[e,,]=x$Err
    Rmse[e,,]=x$Rmse
    Imp[e,]=x$Imp
    Neg[e,]=x$Neg
    Rmax[e,]=x$Rmax
    Psimin[e,]=x$Psimin
    Comm[e,]=x$Comm
    Mu[e,]=x$Mu
    Wald[e,]=x$Wald
    Sec[e,]=x$Sec
  }
  saveRDS(list(design=design,pops=pops,r_pop=r_pop,negpsi_pop=negpsi_pop,mu_pop=mu_pop,structural=structural,
               Hey=Hey,Conv=Conv,Agree=Agree,Err=Err,Rmse=Rmse,Imp=Imp,Neg=Neg,Rmax=Rmax,
               Psimin=Psimin,Comm=Comm,Mu=Mu,Wald=Wald,Sec=Sec),"Simulation_results.rds")

  # Structural populations with their population residuals, uniquenesses and multipliers
  print(pops[structural,])
  print(r_pop[structural])
  print(-negpsi_pop[structural])
  print(mu_pop[structural])

  # Heywood rate of E1-E4 by cell
  print(cbind(design[,c(keys,"n")],apply(Hey,c(1,3),mean)))
  print(sum(Sec)/3600)
}
