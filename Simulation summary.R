# Simulation summary.R
# Tables 2-4 and Figs. 1-3 of Sect. 6.3-6.5 from Simulation_results.rds (Simulation.R).

set.seed(123)
B_auc=1000
res=readRDS("Simulation_results.rds")
design=res$design
nval=sort(unique(design$n))
ncell=nrow(design)

# Groups: sampling (T0, T1, M1, M2, M2s), structural doublets (M2*, M2s*), overshoot (B)
st=res$structural[design$pop]
group=design$truth
group[st & design$truth!="B"]=paste0(design$truth[st & design$truth!="B"],"*")
samp_groups=c("T0","T1","M1","M2","M2s")
str_groups=c("M2*","M2s*","B")
groups=c(samp_groups,str_groups)
hw=res$Hey[,,3] %in% TRUE
dim(hw)=dim(res$Rmax)

# E1, E2 and E3 flag the same samples
agree=res$Hey[,,1]==res$Hey[,,2] & res$Hey[,,2]==res$Hey[,,3]
c(sum(agree,na.rm=TRUE),length(agree))
table(design$m[row(agree)[which(!agree)]],design$n[row(agree)[which(!agree)]])

# Heywood frequency by group and n (Table 2); E4 counts psi <= 0.0051
thr=0.005+1e-4
Freq=matrix(NA,8,5,dimnames=list(groups,nval))
Se=Freq
Ml=Freq
LsThr=Freq
Reps=Freq
for (g in groups)
{
  for (e in 1:5)
  {
    cells=which(group==g & design$n==nval[e])
    x=res$Hey[cells,,3]
    Reps[g,e]=sum(!is.na(x))
    Freq[g,e]=mean(x,na.rm=TRUE)
    Se[g,e]=sqrt(Freq[g,e]*(1-Freq[g,e])/Reps[g,e])
    Ml[g,e]=mean(res$Hey[cells,,4],na.rm=TRUE)
    LsThr[g,e]=mean(res$Psimin[cells,]<=thr,na.rm=TRUE)
  }
}
# Table 2 (ties such as 564/8000 = 0.0705 rounded half up)
floor(Freq*1000+0.5)/1000
max(Se)
round(Ml-LsThr,3)
range(Ml-LsThr)

# Share of T1 among the sampling Heywood cases
for (e in 1:5)
{
  cells=which(!st & design$n==nval[e])
  t1=which(group=="T1" & design$n==nval[e])
  print(c(nval[e],sum(hw[t1,])/sum(hw[cells,])))
}

# Communality change E3 vs E2 and residual among Heywood samples; recovery
MedComm=Freq
MedR=Freq
for (g in groups)
{
  for (e in 1:5)
  {
    cells=which(group==g & design$n==nval[e])
    h=which(hw[cells,,drop=FALSE])
    MedComm[g,e]=median(res$Comm[cells,,drop=FALSE][h])
    MedR[g,e]=median(res$Rmax[cells,,drop=FALSE][h])
  }
}
round(MedComm,3)
round(MedR,3)
cells=which(group=="B" & design$n==1000)
c(E1=median(res$Err[cells,,1]),E2=median(res$Err[cells,,2]),E3=median(res$Err[cells,,3]))

# Fig. 1
cols=c(T0="#0072B2",T1="#56B4E9",M1="#009E73",M2="#E69F00",M2s="#F0E442","M2*"="#D55E00","M2s*"="#CC79A7",B="black")
pchs=c(T0=19,T1=15,M1=17,M2=18,M2s=25,"M2*"=19,"M2s*"=15,B=17)
cairo_pdf("Fig1.pdf",width=7.5,height=3.6)
par(mfrow=c(1,2),mar=c(4,4.5,2,1),mgp=c(2.6,0.7,0),las=1)
for (pn in 1:2)
{
  if (pn==1) gs=samp_groups else gs=str_groups
  plot(NA,xlim=c(50,1000),ylim=c(0,1),log="x",xaxt="n",xlab="Sample size n",
       ylab="Proportion of samples with a Heywood case",
       main=c("(a) Sampling populations","(b) Structural populations")[pn],font.main=1,cex.main=1)
  axis(1,at=nval)
  for (g in gs)
  {
    lines(nval,Ml[g,],lty=3,col=cols[g])
    arrows(nval,Freq[g,]-1.96*Se[g,],nval,Freq[g,]+1.96*Se[g,],angle=90,code=3,length=0.02,col=cols[g])
    lines(nval,Freq[g,],type="b",pch=pchs[g],col=cols[g],bg=cols[g],lwd=1.5)
  }
  legend(c("topright","bottomright")[pn],legend=gs,col=cols[gs],pch=pchs[gs],pt.bg=cols[gs],lwd=1.5,bty="n",cex=0.8)
}
dev.off()

# Tails of the two reports among Heywood samples (Fig. 2)
probs=c(0.5,0.9,0.99,1)
Qr=array(NA,c(2,5,4))
Qe=Qr
for (cl in 1:2)
{
  for (e in 1:5)
  {
    cells=which(st==(cl==2) & design$n==nval[e])
    h=which(hw[cells,,drop=FALSE])
    Qr[cl,e,]=quantile(res$Rmax[cells,,drop=FALSE][h],probs)
    Qe[cl,e,]=quantile(res$Neg[cells,,drop=FALSE][h],probs)
  }
}
max(res$Rmax[hw])
Qr[,,1]
Qe[,,3]
Qe[,,4]
cairo_pdf("Fig2.pdf",width=7.5,height=4.4)
layout(matrix(c(1,2,3,3),2,byrow=TRUE),heights=c(1,0.2))
par(mar=c(4,4.5,2,1),mgp=c(2.6,0.7,0))
for (cl in 1:2)
{
  plot(NA,xlim=c(50,1000),ylim=c(1e-3,1e6),log="xy",xaxt="n",yaxt="n",xlab="Sample size n",
       ylab="Quantile among Heywood samples",
       main=c("(a) Sampling Heywood cases","(b) Structural Heywood cases")[cl],font.main=1,cex.main=1)
  axis(1,at=nval)
  axis(2,at=10^(-3:6),labels=parse(text=paste0("10^",-3:6)),las=1)
  for (q in 1:4)
  {
    lines(nval,pmax(Qr[cl,,q],1e-8),type="b",pch=19,cex=0.6,lty=c(1,2,4,3)[q],col="#D55E00",lwd=1.5)
    lines(nval,pmax(Qe[cl,,q],1e-8),type="b",pch=15,cex=0.6,lty=c(1,2,4,3)[q],col="grey35",lwd=1.5)
  }
}
par(mar=c(0,0,0,0))
plot.new()
legend("center",c("profile residual r_max (E3)","negative uniqueness (E1)","median","90%","99%","max"),
       col=c("#D55E00","grey35",rep("black",4)),lty=c(1,1,1,2,4,3),pch=c(19,15,NA,NA,NA,NA),
       ncol=3,bty="n",cex=0.8)
dev.off()

# Accuracy of each report against its population value (Table 3). mu is the largest multiplier
# of the cap in the capped fit E2; it is zero in a sample without a Heywood case.
tr=matrix(res$r_pop[design$pop],ncell,ncol(res$Rmax))
te=matrix(res$negpsi_pop[design$pop],ncell,ncol(res$Rmax))
tm=matrix(res$mu_pop[design$pop],ncell,ncol(res$Rmax))
er=res$Rmax-tr
ee=res$Neg-te
em=res$Mu-tm
Acc=matrix(NA,5,10,dimnames=list(nval,c("rel_r","rel_mu","rel_e1","half_r","half_mu","half_e1",
                                          "rmse_r","rmse_mu","rmse_r_samp","rmse_mu_samp")))
RmseE1=matrix(NA,5,2,dimnames=list(nval,c("rmse_e1","rmse_e1_samp")))
for (e in 1:5)
{
  c1=which(st & design$n==nval[e])
  c0=which(!st & design$n==nval[e])
  Acc[e,1]=median(abs(er[c1,])/tr[c1,])
  Acc[e,2]=median(abs(em[c1,])/tm[c1,])
  Acc[e,3]=median(abs(ee[c1,])/te[c1,])
  Acc[e,4]=mean(abs(er[c1,])<=0.5*tr[c1,])
  Acc[e,5]=mean(abs(em[c1,])<=0.5*tm[c1,])
  Acc[e,6]=mean(abs(ee[c1,])<=0.5*te[c1,])
  Acc[e,7]=sqrt(mean(er[c1,]^2))
  Acc[e,8]=sqrt(mean(em[c1,]^2))
  Acc[e,9]=sqrt(mean(er[c0,]^2))
  Acc[e,10]=sqrt(mean(em[c0,]^2))
  RmseE1[e,1]=sqrt(mean(ee[c1,]^2))
  RmseE1[e,2]=sqrt(mean(ee[c0,]^2))
}
signif(Acc,3)
# RMSE of the E1 uniqueness, quoted in the text only (its population value is infinite for p = 3)
signif(RmseE1,3)
range(res$mu_pop[res$structural])

# Separation of structural from sampling Heywood cases: pairwise AUCs (Table 4, Fig. 3).
# Four reports: residual r_max (E3), multiplier mu (E2), negative uniqueness (E1), and the Wald
# statistic of the smallest unrestricted uniqueness with its sign changed. A Heywood sample
# without a Wald value (no convergence) counts as zero, the value least in favor of "structural".
auc=function(score,y)
{
  n1=sum(y)
  n0=sum(!y)
  return((sum(rank(score)[y])-n1*(n1+1)/2)/(n1*n0))
}

wz=-res$Wald
wz[is.na(wz)]=0
scores=list(r=res$Rmax,mu=res$Mu,e1=res$Neg,wald=wz)
nsc=length(scores)
# Share of Heywood samples without a Wald value, by n
for (e in 1:5)
{
  cells=which(design$n==nval[e])
  print(c(nval[e],mean(is.na(res$Wald[cells,,drop=FALSE][hw[cells,,drop=FALSE]]))))
}

pairs=expand.grid(s0=samp_groups,s1=str_groups,stringsAsFactors=FALSE)[,2:1]
npair=nrow(pairs)
Auc=array(NA,c(nsc,npair,5),dimnames=list(names(scores),paste(pairs$s1,"vs",pairs$s0),nval))
Sc=list()
for (e in 1:5)
{
  Sc[[e]]=list()
  for (g in groups)
  {
    cells=which(group==g & design$n==nval[e])
    h=which(hw[cells,,drop=FALSE])
    Sc[[e]][[g]]=matrix(NA,length(h),nsc)
    for (s in 1:nsc)
    {
      Sc[[e]][[g]][,s]=scores[[s]][cells,,drop=FALSE][h]
    }
  }
  for (j in 1:npair)
  {
    a=pairs$s1[j]
    b=pairs$s0[j]
    if (nrow(Sc[[e]][[a]])>=20 && nrow(Sc[[e]][[b]])>=20)
    {
      y=rep(c(TRUE,FALSE),c(nrow(Sc[[e]][[a]]),nrow(Sc[[e]][[b]])))
      for (s in 1:nsc)
      {
        Auc[s,j,e]=auc(c(Sc[[e]][[a]][,s],Sc[[e]][[b]][,s]),y)
      }
    }
  }
}
Ar=Auc["r",,]
Ae=Auc["e1",,]
D=Ar-Ae
round(D,2)
round(Ar-Auc["mu",,],2)
round(Ar-Auc["wald",,],2)

# Heywood samples per group and n
Cnt=matrix(NA,8,5,dimnames=list(groups,nval))
for (e in 1:5)
{
  for (g in groups)
  {
    Cnt[g,e]=nrow(Sc[[e]][[g]])
  }
}
Cnt

# Pairs with at least 50 Heywood samples in each group at every n (Table 4), their average and
# stratified bootstrap intervals. The same resamples are used for the four reports, so the
# intervals for the differences are paired.
common=which(apply(Cnt[pairs$s1,],1,min)>=50 & apply(Cnt[pairs$s0,],1,min)>=50)
rownames(Ar)[common]
# For comparison: average over the pairs with at least 20 Heywood samples at every n
common20=which(rowSums(is.na(Ar))==0)
rownames(Ar)[common20]
round(apply(Auc[,common20,],c(1,3),mean),3)
rows=c("r","mu","e1","wald","r-mu","r-e1","r-wald")
Avg=matrix(NA,21,5,dimnames=list(paste0(rep(rows,each=3),c("","_lo","_hi")),nval))
for (e in 1:5)
{
  bs=matrix(NA,B_auc,nsc)
  for (b in 1:B_auc)
  {
    rs=list()
    for (g in groups)
    {
      i=sample.int(nrow(Sc[[e]][[g]]),nrow(Sc[[e]][[g]]),replace=TRUE)
      rs[[g]]=Sc[[e]][[g]][i,,drop=FALSE]
    }
    x=matrix(NA,length(common),nsc)
    for (jj in 1:length(common))
    {
      a=pairs$s1[common[jj]]
      g0=pairs$s0[common[jj]]
      y=rep(c(TRUE,FALSE),c(nrow(rs[[a]]),nrow(rs[[g0]])))
      for (s in 1:nsc)
      {
        x[jj,s]=auc(c(rs[[a]][,s],rs[[g0]][,s]),y)
      }
    }
    bs[b,]=colMeans(x)
  }
  est=rowMeans(Auc[,common,e])
  col=vector()
  for (s in 1:nsc)
  {
    col=c(col,est[s],quantile(bs[,s],c(0.025,0.975)))
  }
  for (s in 2:nsc)
  {
    col=c(col,est[1]-est[s],quantile(bs[,1]-bs[,s],c(0.025,0.975)))
  }
  Avg[,e]=col
}
round(Avg,3)

# Fig. 3
pal=colorRampPalette(c("#053061","#2166AC","#92C5DE","#F7F7F7","#F4A582","#B2182B","#67001F"))(101)
lim=max(abs(D),na.rm=TRUE)
cairo_pdf("Fig3.pdf",width=8,height=4)
layout(matrix(1:3,1,3),widths=c(1,1.15,0.12))
par(mar=c(4,4.5,2,1))
plot(NA,xlim=c(50,1000),ylim=c(0.45,1),log="x",xaxt="n",xlab="Sample size n",ylab="AUC",
     main="(a) Pairs with 50 or more cases at every n",font.main=1,cex.main=1)
axis(1,at=nval)
abline(h=0.5,lty=3,col="grey50")
for (j in common)
{
  lines(nval,Ar[j,],col=adjustcolor("#D55E00",0.4))
  lines(nval,Ae[j,],col=adjustcolor("grey35",0.4))
}
lines(nval,Avg["mu",],type="b",pch=17,col="#0072B2",lwd=2.5)
lines(nval,Avg["wald",],type="b",pch=18,col="#009E73",lwd=2.5)
lines(nval,Avg["e1",],type="b",pch=15,col="grey35",lwd=2.5)
lines(nval,Avg["r",],type="b",pch=19,col="#D55E00",lwd=2.5)
legend("bottomright",c("r_max (E3)","multiplier (E2)","negative uniqueness (E1)","Wald statistic"),
       col=c("#D55E00","#0072B2","grey35","#009E73"),pch=c(19,17,15,18),lwd=2.5,bty="n",cex=0.85)
par(mar=c(4,7,2,0.5))
plot(NA,xlim=c(0.5,5.5),ylim=c(0.5,npair+0.5),axes=FALSE,xlab="Sample size n",ylab="",
     main="(b) AUC difference, r_max minus E1",font.main=1,cex.main=1)
for (j in 1:npair)
{
  for (e in 1:5)
  {
    yj=npair+1-j
    if (is.na(D[j,e]))
    {
      rect(e-0.5,yj-0.5,e+0.5,yj+0.5,col="grey92",border=NA)
    } else {
      ci=round((D[j,e]+lim)/(2*lim)*100)+1
      rect(e-0.5,yj-0.5,e+0.5,yj+0.5,col=pal[ci],border=NA)
      text(e,yj,sprintf("%+.2f",D[j,e]),cex=0.7,col=ifelse(abs(D[j,e])>0.6*lim,"white","black"))
    }
  }
}
axis(1,at=1:5,labels=nval,tick=FALSE)
axis(2,at=npair:1,labels=rownames(D),las=1,tick=FALSE,cex.axis=0.8)
par(mar=c(4,0.5,2,2.5))
image(1,seq(-lim,lim,length.out=101),matrix(seq(-lim,lim,length.out=101),1),col=pal,axes=FALSE,xlab="",ylab="")
axis(4,las=1,cex.axis=0.8)
box()
dev.off()
