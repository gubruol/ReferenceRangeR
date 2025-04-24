##################################################################
# TML - Reference limit estimator                              ###
# Author: F. Arzideh                                           ###
# https://www.degruyter.com/document/doi/10.1515/CCLM.2007.249/html
# 08.05.2024 Slight modifications by G. Brandhorst             ###
##################################################################

### Modifications from original functions ###
# gr_est0.R: line 31: return(list(P=res$P,mu=res$mu,sig=res$sig,T2=res$T2,DL975=res$DL975,DL95=res$DL95,
# myplot=res$myplot,DL25=res$DL25))
# gr_all0.R: line 135: myplot <- recordPlot()
# gr_all0.R: line 136: return(list(P=p,mu=mu,sig=sig,T2=T2,DL975=DL975,DL95=DL950,myplot=myplot,DL25=DL25))
# legend removed: legend(DL950 ... 2 lines
# mtext removed: mtext(paste("Estimated distributions ... 3 lines
# removed: #mtext(date(),cex=0.6,font=3,outer=F,line=5,side=1,adj=1)

##################################################################
# DL_pnm function calculates RLs for power-normal distribution ###
# updated: 29.09.2014 (without DL)                             ###
# Autor: Arzideh                                               ###
# 26.06.2019 (Arzideh): for different percentiles              ###
##################################################################
#library(msm)
#library(geoR)
#####################################
DL.pnm<-
function (mu0,sig0,p,lam,TL,pc1=0.025,pc2=0.05,PC1=1-pc1,PC2=1-pc2) 
{
# DL.pnm function: 22.04.09
# p, mu0, sig0 and lam are the estimated parameters for the distr. of non-pathological values   
########################################
dl975<-qtnorm(PC1,mu0,sig0,lower=TL)
dl950<-qtnorm(PC2,mu0,sig0,lower=TL)
density975<-dtnorm(dl975,mu0,sig0,lower=TL)
density950<-dtnorm(dl950,mu0,sig0,lower=TL)
DL975<-(dl975*lam+1)^(1/lam)
DL950<-(dl950*lam+1)^(1/lam)
Den975<-(density975)*(DL975^(lam-1))
Den950<-(density950)*(DL950^(lam-1))
####################################
dl25<-qtnorm(pc1,mu0,sig0,lower=TL)
dl5<-qtnorm(pc2,mu0,sig0,lower=TL)
density25<-dtnorm(dl25,mu0,sig0,lower=TL)
density5<-dtnorm(dl5,mu0,sig0,lower=TL)
DL25<-(dl25*lam+1)^(1/lam)
DL5<-(dl5*lam+1)^(1/lam)
Den25<-(density25)*(DL25^(lam-1))
Den5<-(density5)*(DL5^(lam-1))
#########################
return(list(DL25=DL25,DL5=DL5,DL975=DL975,DL950=DL950,
Den25=Den25,Den5=Den5,Den975=Den975,Den950=Den950,mu=mu0,sigma=sig0,P=p))
}


################################################################
# DL_log function calculates RLs for log-normal distribution ###
# updated: 29.09.2014 (without DL)                           ###
# Autor: Arzideh                                             ###
# 26.06.2019 (Arzideh): with different percentiles           ### 
################################################################
#library(msm)
#library(geoR)
##############
DL.log<-
function (mu0,sig0,p,pc1=0.025,pc2=0.05,PC1=1-pc1,PC2=1-pc2) 
{
# DL.log function: 25.01.10
# p, mu0, sig0 and are the estimated parameters for the distr. of non-pathological values   
########################################
dl975<-qnorm(PC1,mu0,sig0)
dl950<-qnorm(PC2,mu0,sig0)
density975<-dnorm(dl975,mu0,sig0)
density950<-dnorm(dl950,mu0,sig0)
DL975<-exp(dl975)
DL950<-exp(dl950)
Den975<-(density975)/(DL975)
Den950<-(density950)/(DL950)
dl25<-qnorm(pc1,mu0,sig0)
dl5<-qnorm(pc2,mu0,sig0)
density25<-dnorm(dl25,mu0,sig0)
density5<-dnorm(dl5,mu0,sig0)
DL25<-exp(dl25)
DL5<-exp(dl5)
Den25<-(density25)/(DL25)
Den5<-(density5)/(DL5)
return(list(DL25=DL25,DL5=DL5,DL975=DL975,DL950=DL950,
Den25=Den25,Den5=Den5,Den975=Den975,Den950=Den950,mu=mu0,sigma=sig0,P=p))
}


# last modification: 07.03.2016 (Farhad Arzideh), if (wv >=1 || is.na(wv))
# est_dtpn1.R
###################################
est.dtpn1<-
function (d=data,l0=0.5,sl=0.1,n=10,TR=max(data),mod=0.8,com=2,c2=1,TL2=round(min(mod-c2*(10^(-com)),mod-(mod-quantile(d,0.05))*0.1),com),
TL1=min(d),I=100,I1=25,T=max(d),nen=1,epsilon1=1E-4,steptol=5) 
{
#################################
# function est.dtpn1 (parameter estimation of TPN distribution to the data truncated at different points) written on 12.08.2010
# for truncated only at the low values (L)
# PN(lambda,mu,sigma)
# l0, sl and n: define the range in which the transformation parameter lambda is searched:
# e.g.: l0=0.5,sl=0.1 and n=10 means proper lambda is searched in {0.6,0.7,0.8,...,1.5}
# TR: fix truncation point at the right side of mod
# com: accuracy of data values; e.g. for data values 2.11, 1.67 is com=2 and for 144,150 is com=0
#  truncation points are determined at the left side of mod of the data in interval: TL=[TL1,TL2]
# c2 defines how far from mod TL2 is searched
# e.g.for creatinin with data given with two decimals(0.94 or 1.03)and mod=0.88 and for c2=10, TL1 and TL2 can be calculated as:
# TL1=min(data),
# TL2= 0.88-10*(10^-2)=0.78, and thereby:
# TL=[TL1,0.78]
# T: limit of data; data >=T are eliminated. Default: T=max(d)
# I and epsilon: fixed variables to stop the numerical procedure if predefined criterion are satisfied.
# nen: correction factor to modify the proper choice for start value to estimate sigma. default: nen=1 
# mod: mod (mod of data), start value for estimation of mu. 
#################################
d<-d[0<d & d<=T] # eliminate extreme values
q2<-(length(d[d<=TR]))/length(d) # q2: fraction of data smaller than TR, which are not used to estimate the parameters
dt<-T1point(data=d,TL1=TL1,TL2=TL2) # delivers truncation points (TL) and fractions of data <=TL (q1).
q1<-dt$q1
T1<-dt$T1
m<-length(q1) # m: number of truncation points, which gives us m different estimations! 
y<-matrix(rep(0,m*13),ncol=13)
y0<-matrix(rep(0,13),ncol=13)
# now for each truncation point from m truncation points the fitted distribution is estimated.
i<-1;
while(i<=m)
         {
         qu1<-q1[i]
         tr1<-T1[i]
         est<-est.tpn(d=d,mu=mod,q2=q2,T2=TR,q1=qu1,T1=tr1,I=I,I1=I1,l0=l0,sl=sl,n=n,
         nen=nen,com=com,epsilon1=epsilon1,steptol=steptol)$h1
         if (est[11]==-1E+7)
         {
         y[i,1:12]<-c(rep(1,12))
         y[i,13]<-1000
         }
         else
        {
         T11<-est[3]
         T21<-est[4]
         p<-est[9]
         mu<-est[6]
         sig<-est[7]
         frac<-(q2-qu1)/est[8]
         lam<-est[10]
         ll<-est[11]
         d1<-d[d>=tr1 & d<=TR]
         if (lam==0)
            {
             d2<-log(d1)
             d3<-log(d)
            crit<-log(TL2)
            }
         else
            {
            d2<-(d1^lam-1)/lam
            d3<-(d^lam-1)/lam
            crit<-((TL2)^lam-1)/lam
            }
         w<-KD.tnorm(data=d2,mu=mu,sig=sig,T1=T11,T2=T21)$K 
         v<-KD.tnorm2(data=d3,p=frac,mu=mu,sig=sig,crit=crit)
         wv<-w+v
         if (wv >=1 || is.na(wv))
         {
         y[i,1:12]<-c(rep(1,12))
         y[i,13]<-0
         } else
        {
         wvinv<-1/wv
         y[i,1:13]<-c(qu1,q2,tr1,TR,p,mu,sig,lam,ll,w,v,v+w,wvinv)
         }
        } 
         i<-i+1;
       }
index0<-which.min(y[,12])
y0[1,]<-y[index0,]
y<-round(y,5)
y0<-round(y0,5)
colnames(y)<-c("q1","q2","T1","T2","P","mu","sigma","lam","LL","KD1","KD2","KD","invKD")
colnames(y0)<-c("q1","q2","T1","T2","P","mu","sigma","lam","LL","KD1","KD2","KD","invKD")
return(list(y=y,y0=y0))
}



# last modification: 07.03.2016 (Farhad Arzideh), if (wv >=1 || is.na(wv))
# est_dtpn2.R
#######################################
est.dtpn2<-
function (d=data,l0=0.5,sl=0.1,n=10,TL=0,mod=0.8,com=2,c2=1,TR2=max(d),TR1=round(max(mod+c2*(10^(-com)),mod+(min(quantile(d,0.95),TR2)-mod)*0.1),com),I=100,I1=25,T=max(d),nen=1,epsilon1=1E-4,steptol=5) 
{
###############################
# last modification: 05.03.2012: (for the case that no MLE delivered: y[,13]<-0).
# function est.dtpn2 (parameter estimation of TPN distribution to the data truncated at different points)  written on 12.08.2010
# for truncated only at high values (H)
# PN(lambda,mu,sigma)
# l0, sl and n: define the range in which the transformation parameter lambda is searched:
# e.g.: l0=0.5,sl=0.1 and n=10 means proper lambda is searched in {0.6,0.7,0.8,...,1.5}
# TL: fix truncation point at the left side of mod
# com: accuracy of data values; e.g. for data values 2.11, 1.67 is com=2 and for 144,150 is com=0
#  truncation points are determined at the right side of mod of the data in interval: TR=[TR1,TR2]
# c2 defines how far from mod TR1 is searched
# e.g. for creatinin with data given with two decimals (0.94 or 1.03) and mod=0.88 and for c2=10, TR1 and TR2 can be calculated as:
# , TR2=max(data),
#  TR1= 0.88+10*(10^-2)=0.98, and thereby:
#  TR=[0.98,TR2]
# T: limit of data; data >=T are eliminated. Default: T=max(d)
# I and epsilon: fixed variables to stop the numerical procedure if predefined criterion are satisfied.
# nen: correction factor to modify the proper choice for start value to estimate sigma. default: nen=1 
# mod: mod (mod of data), start value for estimation of mu. 
##############################
d<-d[d<=T]                                                  # eliminate extreme values
q1<-(length(d[d<TL]))/length(d)                             # q1: fraction of data smaller than TL, which are not used to estimate the parameters
dt<-T2point(data=d,TR1=TR1,TR2=TR2)              # delivers truncation points (T2) and fractions of data <=T2 (q2).
q2<-dt$q2
T2<-dt$T2
m<-length(q2)                                      # m: number of truncation points, which gives us m different estimations! 
y<-matrix(rep(0,m*13),ncol=13)
y0<-matrix(rep(0,13),ncol=13)
################################
# now for each truncation point from m truncation points the fitted distribution is estimated.
i<-1;
while(i<=m)
       {
         qu2<-q2[i]
         tr2<-T2[i]
         est<-est.tpn(d=d,mu=mod,q1=q1,T1=TL,q2=qu2,T2=tr2,I=I,I1=I1,l0=l0,sl=sl,n=n,
         nen=nen,com=com,epsilon1=epsilon1,steptol=steptol)$h1
         if (est[11]==-1E+7)
         {
         y[i,1:12]<-c(rep(1,12))
         y[i,13]<-0
         }
         else
        {
         T11<-est[3]
         T21<-est[4]
         p<-est[9]
         mu<-est[6]
         sig<-est[7]
         frac<-(qu2-q1)/est[8]
         lam<-est[10]
         ll<-est[11]
         d1<-d[d>=TL & d<=tr2]
         if (lam==0)
            {
             d2<-log(d1)
             d3<-log(d) 
            crit<-log(TR1)
            }
         else
            {
            d2<-(d1^lam-1)/lam
            d3<-(d^lam-1)/lam
            crit<-((TR1)^lam-1)/lam
            }
         w<-KD.tnorm(data=d2,mu=mu,sig=sig,T1=T11,T2=T21)$K 
         v<-KD.tnorm1(data=d3,p=frac,mu=mu,sig=sig,T1=T11,crit=crit)
         wv<-w+v
         if (wv >=1 || is.na(wv))
         {
         y[i,1:12]<-c(rep(1,12))
         y[i,13]<-0
         } else
        {
         wvinv<-1/wv
         y[i,1:13]<-c(q1,qu2,TL,tr2,p,mu,sig,lam,ll,w,v,v+w,wvinv)
         }
        } 
         i<-i+1;
       }
index0<-which.min(y[,12])
y0[1,]<-y[index0,]
y<-round(y,5)
y0<-round(y0,5)
colnames(y)<-c("q1","q2","T1","T2","P","mu","sigma","lam","LL","KD1","KD2","KD","invKD")
colnames(y0)<-c("q1","q2","T1","T2","P","mu","sigma","lam","LL","KD1","KD2","KD","invKD")
return(list(y=y,y0=y0))
}

EST.PN.L<-function(data,Q1,Q2,dec=dec.num,mod,l0=0,n=n,dif1,nen=1,lstep=0.1,nstep=21)
 {   
# F. Arzideh
# 15.08.2019
# estimate of parameter of power normal distribution
# this function returns 'est': estimated parameters and 'estim': controls whether parameters have been estimated:
# estim=0 means no estimation, and estim=1 means estimation.
# two other required var. for the next calculation (for CI) are also returned: q1 and maxT 
##################################
##################################
TL1<-round(quantile(data,Q1),dec)
TR<-round(quantile(data,Q2),dec) 
est<-est.dtpn1(d=data,l0=l0,sl=0.1,n=n,TR=TR,mod=mod,com=dec,c2=dif1,TL1=TL1,nen=nen)
est1<-est$y0
#########################
   qstep<-0.01
   q20<-Q2
   while (est1[1,12]==1 && qstep <= 1-q20) {
      q20<-round(q20+qstep,2)
      TR<-round(quantile(data,q20),dec)
      est<-est.dtpn1(d=data,l0=l0,sl=0.1,n=n,TR=TR,mod=mod,com=dec,c2=dif1,TL1=TL1,nen=nen)
      est1<-est$y0
      }
################################
   if (est1[1,12]==1) 
            {# if no proper truncation points are found the program is stopped!!
             estim<-0
             }else{
             estim<-1
             TL<-est1[3]
             TR<-est1[4]
             ## final estimation (accuaracy of lambda: 0.01):
             est<-est.dtpn1(d=data,l0=max(0,est1[8]-lstep),sl=0.01,n=nstep,TR=TR,mod=mod,
             com=dec,TL1=TL,TL2=TL,nen=nen)
                }
return(list(est=est,estim=estim,q2=q20))
}

EST.PN<-function(data,Q1,Q2,dec=dec.num,mod,l0=0,n=n,dif1,nen=1,lstep=0.1,nstep=21)
 {   
# F. Arzideh
# 13.08.2019
# estimate of parameter of power normal distribution
# this function returns 'est': estimated parameters and 'estim': controls whether parameters have been estimated:
# estim=0 means no estimation, and estim=1 means estimation.
# two other required var. for the next calculation (for CI) are also returned: q1 and maxT 
# 10.12.2019: max.TR > mod was modified as max.TR > med.emp (median of data) and output maxT as calcvulated TR2. 
##################################
##################################
max.TR<-4*(mod-min(data))+mod
med.emp<-quantile(data,0.5)
TL<-round(quantile(data,Q1),dec)
if (max.TR > med.emp)
     {
      TR2<-round(min(quantile(data,Q2),max.TR),dec)
      }else{
      TR2<-round(quantile(data,Q2),dec) 
      } 
est<-est.dtpn2(d=data,l0=l0,sl=0.1,n=n,TL=TL,mod=mod,com=dec,c2=dif1,TR2=TR2,nen=nen)
est1<-est$y0
#########################
   qstep<-0.01
   q10<-Q1
   while (est1[1,12]==1 && qstep <=q10) {     # if the truncated part is too narrow, it will be boarder.
      q10<-round(q10-qstep,2)
      TL<-round(quantile(data,q10),dec)
      est<-est.dtpn2(d=data,l0=l0,sl=0.1,n=n,TL=TL,mod=mod,com=dec,c2=dif1,TR2=TR2,nen=nen)
      est1<-est$y0
      }
################################
   if (est1[1,12]==1) 
            {# if no proper truncation points are found the program is stopped!!
             estim<-0
             }else{
             estim<-1
             TR<-est1[4]
             TL<-est1[3]
             ## final estimation (accuaracy of lambda: 0.01):
             est<-est.dtpn2(d=data,l0=max(0,est1[8]-lstep),sl=0.01,n=nstep,TL=TL,mod=mod,
             com=dec,TR1=TR,TR2=TR,nen=nen)
                }
return(list(est=est,estim=estim,q1=q10,maxT=TR2))
}


#source("TNorm.txt")
##################

##################
est.tpn<-function (d,q1=0.0,q2,T1=quantile(d,q1),T2=quantile(d,q2),com,mu=128,I1=50,I=200,l0,sl,n,
epsilon1=1E-4,nen=1,steptol=4)
 {
####################################
  # last modification: 30.05.2012: (h2 eliminated.) 
  # this function evaluate the fitted truncated Power-Normal distribution
  # modified on 05.09.2019: start values are modified:
  # for log-ND: mu_Y = ln(mod_X)+ (sigma_Y)**2
  # for PN:     mu_Y = BC(mod_X, lambda) + {(1-lambda)*(sigma_Y)**2}/ (mod_X)**lambda
  # where X non-transformed, Y transformed, mu_Y and sigma_Y= mu and sigma of ND,
  # mod_X: mode of non-transformed vlaues are.     
###################################
  er1<-0.5*(10^(-com)) 
  T1<-max(er1,T1-er1)
  T2<-T2+er1
  dT<-d[d>T1 & d<T2]
  h<-matrix(rep(0,n*14),nrow=n)
  h1<-matrix(rep(0,14),nrow=1)
  i1<-0; while(i1<=n-1)
  {
  lam<-l0+i1*sl;
  if (lam==0)
            {
            T11<-log(T1)
            T21<-log(T2)
            d1<-log(d)
            d2<-log(dT)
            sig<-nen*(sd(d2))/(q2)
            #mu1<-log(mu)
            mu1<-log(mu)+sig^2
            } else
       {
       T11<-(T1^lam-1)/lam 
       T21<-(T2^lam-1)/lam
       d1<-(d^lam-1)/lam
       d2<-(dT^lam-1)/lam
       sig<-nen*(sd(d2))/(q2)
       mu1<-(mu^lam-1)/lam + (1-lam)*(sig^2)/(mu^lam)
       #mu1<-(mu^lam-1)/lam 
       }
####################################
  u<-(lam-1)*(sum(log(dT)))
  #sig<-nen*(sd(d2))/(q2)
  v<-TNorm(q1=q1,q2=q2,T1=T11,T2=T21,d=d1,mu=mu1,sig=sig,I=I1,er1=0,epsilon1=1E-2)
  step<-1
  while (v$differential > 0.01 && (step < steptol))
        {              
         nen<-1-(step/steptol)
         sig<-nen*(sd(d2))/(q2)
         v<-TNorm(q1=q1,q2=q2,T1=T11,T2=T21,d=d1,mu=mu1,sig=sig,I=I1,er1=0,epsilon1=1E-2)
         step<-step+1        
         }
#######################################
  step0<-1
  while (v$differential > 0.01 && (step0 < steptol))
        {              
         nen<-1+(step0/steptol)
         sig<-nen*(sd(d2))/(q2)
         v<-TNorm(q1=q1,q2=q2,T1=T11,T2=T21,d=d1,mu=mu1,sig=sig,I=I1,er1=0,epsilon1=1E-2)
         step0<-step0+1         
         }
######################################
  i<-i1+1
  if (v$differential > 0.01)
                {
                 h[i,1]<-h[i,2]<-h[i,3]<-h[i,4]<-h[i,5]<-h[i,6]<-h[i,7]<-h[i,8]<-h[i,9]<-h[i,10]<-h[i,12]<-h[i,13]<-1E+7
                 h[i,11]<- -1E+7
                 h[i,14]<- 1E+7
                 } else
                {

  v<-TNorm(q1=q1,q2=q2,T1=T11,T2=T21,d=d1,mu=mu1,sig=sig,I=I,er1=0,epsilon1=epsilon1)
  if (v$differential > epsilon1)
                {
                 h[i,1]<-h[i,2]<-h[i,3]<-h[i,4]<-h[i,5]<-h[i,6]<-h[i,7]<-h[i,8]<-h[i,9]<-h[i,10]<-h[i,12]<-h[i,13]<-1E+7
                 h[i,11]<- -1E+7
                 h[i,14]<- 1E+7
                 } else
                {
  h[i,1]<-round(v$q1,8)
  h[i,2]<-round(v$q2,8)
  h[i,3]<-round(v$T1,8)
  h[i,4]<-round(v$T2,8)
  h[i,5]<-round(v$I,4)
  h[i,6]<-round(v$mu,8)
  h[i,7]<-round(v$sig,8)
  h[i,8]<-round(v$Phi,8)
  h[i,9]<-round(v$frac,4)
  h[i,10]<-lam
  h[i,11]<-round(v$l+u,4)
  h[i,12]<-round(T2,4)
  h[i,13]<-round(T1,4)
  h[i,14]<-round(v$differential,4)
               }}
  i1<-i1+1
  }
  index<-which.max(h[,11])
  h1<-h[index,]
  return(list(h=h,h1=h1))
}


##################################################################
# - gr.all0 function displays the estimated distributions for  ###
#   non-pathological and pathological values for path.=B       ###
# - updated: 29.09.2014 (without DL) and some lables           ###
# - updated: 08.10.2014 legends and footnotes are modified     ###
#   (Autor: Arzideh)                                           ###
# - 02.02.2015 (Ar): modified to:                              ###
#    i) display path. values on both sides (low and high       ### 
#       values),                                               ###
#    ii) display the estimated distributions for all 3         ###
#        cases (H, B, L),                                      ###
#        functions gr_all1 and gr_all2 have been eliminated    ###
#    iii) apply the tool as default in windows-system          ###
#         or alternatively in linux-system                     ###
# - 05.02.2015 (Ar): 5% and 95% were removed                   ###
# - 27.07.2016 (Arzideh)  y-limits has been fitted             ###
# - 17.08.2016 (Arzideh): for graphical purpose :              ###
#   kden<-kde1(data=d,bw=bwd,nb=nb,mod=mod,q1=0,x3=X3)         ###
# - 26.06.2019 (Arzideh): different percentiles:               ###
#                          pc1, pc2,PC1,PC2                    ###
##################################################################
#library(msm) 
#library(geoR) 
#source("kde.txt")
#source("kde1.txt")
#source("DL_pnm.txt")
################################
################################
gr.all0<-
function (d,lam=0.15,q1=0.025,T1=quantile(d,q1),q2=0.64,T2=quantile(d,q2),sig=34.5,mu=129.28,p=0.93,lam2=lam2,
eps1=0.2,x12=0.01,over=30000,er1,x01=x01,x02=x02,step0=bw,nb=2048,c1=0.8,c2=1,sc1=1.2,
labx=labx,fac1=0.75,bwd=bw,low,high,com=2,main="estimation",pc1=0.025,pc2=0.05,PC1=0.975,PC2=0.95) 
{
# written at 05th March 2009
# last modification: 22.08.2012, border: gray15
# for both-sides truncated (B)
# gr.all0 evaluates the DLs and displays the estimated distributions
#########################################
kdew<-kde(data=d,nb=nb,q1=0,bw0=bwd)
y00<-kdew$k
mod<-kdew$u
y00<-round(y00*sc1,3)
minD<-min(d)
maxD<-max(d)
if(maxD > 2*minD) X3<-0 else X3<-minD 
kden<-kde1(data=d,bw=bwd,nb=nb,mod=mod,q1=0,x3=X3)
# kde calculates a kernel density estimation with Sheather & Jones Method 
k1<-kden$k1
u1<-kden$u1 
# k1 and u1 give kde in all points
k2<-kden$k2
u2<-kden$u2
m2<-kden$m2
# k2 and u2 give kde in all points >=mode and m2 length of k2 
k3<-kden$k3
u3<-kden$u3
m3<-kden$m3
# k3 and u3 give kde in all points < mode and m3 length of k3 
bw1<-kden$bw1
k4<-p*(dboxcox(u2,lambda=lam,mean=mu,sd=sig)) # density of pn estimated for >=mod 
te<-rep(0,m2)
j<-1
while (j<=m2)
         {
          te[j]<-max(k2[j]-k4[j],0)
          j<-j+1
          }
te<-te
fx1<-p*(dboxcox(u1,lambda=lam,mean=mu,sd=sig))
###########################################
if (lam==0)
                 {
                 dl<-DL.log(mu0=mu,sig0=sig,p=p,pc1=pc1,pc2=pc2,PC1=PC1,PC2=PC2)
                 }
else
                 {
                dl<-DL.pnm(mu0=mu,sig0=sig,p=p,lam=lam,TL=-1/lam,pc1=pc1,pc2=pc2,PC1=PC1,PC2=PC2)
                 }
############################################
k5<-p*(dboxcox(u3,lambda=lam,mean=mu,sd=sig)) # density of pn estimated for < mod 
tet<-rep(0,m3)
j<-1
while (j<=m3)
         {
          tet[j]<-max(k3[j]-k5[j],0)
          j<-j+1
          }
tet<-tet
#####################################
DL25<-round(dl$DL25,com)+lam2
DL5<-round(dl$DL5,com)+lam2
DL975<-round(dl$DL975,com)+lam2
DL950<-round(dl$DL950,com)+lam2
Den25<-dl$Den25
Den95<-dl$Den950
Den975<-dl$Den975
Den5<-dl$Den5
f<-round(seq(0,y00,y00/4),3)
s<-(step0)*3
f1<-seq(x01,x02,s)
d<-d+lam2
u1<-u1+lam2
u2<-u2+lam2
u3<-u3+lam2
density02<-u2[te==0 & u2<=high]
P.min<-density02[which.min(high-density02)]
density03<-u3[tet==0 & u3>=low]
P.max<-density03[which.min(density03-low)]
# windows(width=7,height=5)
op<-par(mgp = c(1, 0.4, 0),mar=c(6,2,2,3))
Y00<-max(hist(d,br=c(-10,seq(x01+0.5*step0,x02+0.5*step0,step0),over),prob=TRUE)$density,fx1)
hist(d,br=c(-10,seq(x01+0.5*step0,x02+0.5*step0,step0),over),prob=TRUE,main=main,xlab=labx,lwd=1,
xlim=c(x01,x02),ylim=c(0,max(y00,Y00)),cex.axis=0.7,cex.lab=0.8,cex.main=0.9,border="gray50",axes=FALSE)
#mtext(paste("Estimated distributions for non-pathological values (green curve), pathological values (red)","\n",
#            "and whole data (blue). Green lines (and given numbers) indicate 2.5 and 97.5 percentiles of","\n",
#             "the estimated distribution for non-pathological values (RL)."),cex=0.8,font=3,outer=F,line=4,side=1)
par(adj = 1)
#mtext(date(),cex=0.6,font=3,outer=F,line=5,side=1,adj=1)
lines(u1,fx1,lwd=2,lty=1,col="green4")
lines(u3[u3<=P.max],tet[u3<=P.max],lwd=0.4,lty=1,col="red")
lines(u2[u2>=P.min],te[u2>=P.min],lwd=0.4,lty=1,col="red")
axis(side=2,at=f,pos=c(x01,0),cex.axis=0.7,lwd=1,tcl=-0.25)
axis(side=1,at=f1,pos=c(0,0),cex.axis=0.7,lwd=1,tcl=-0.25)
lines(density(x=d,kernel="gaussian",n=nb,cut=3,from=x01,bw=bwd),col="blue",lty=2)
lines(x=c(DL975,DL975),y=c(0,Den975),col="green4",lty=1,lwd=2)
lines(x=c(DL25,DL25),y=c(0,Den25),col="green4",lty=1,lwd=2)
text(DL25,3*Den25,font=1,col="green4",cex=0.8,pos=2,labels=eval(DL25))
text(DL975,3*Den975,font=1,col="green4",cex=0.8,pos=4,labels=eval(DL975))
##################
#legend(DL950,c2*0.95*y00,legend=c("whole data set","non-pathological","pathological"),col=c("blue","green4","red"),
#lty=c(2,1,1),cex=0.7,lwd=1,bty="n")
abline(0,0)
par(adj =0.5)
par(op)
myplot <- recordPlot()
return(list(P=p,mu=mu,sig=sig,T2=T2,DL975=DL975,DL95=DL950,myplot=myplot,DL25=DL25))
}


################################################################
# gr.est0 function displays the estimated distributions for  ###
# non-pathological and pathological values for path.=B       ###
# updated: 29.09.2014 (without DL) and some lables           ###
# Autor: Arzideh                                             ###
# 02.02.2015 (Ar): modified to:                              ###
#    i) display path. values on both sides (low and high     ### 
#       values),                                             ###
#    ii) display the estimated distributions for all 3       ###
#        cases (H, B, L),                                    ###
#        functions gr_est1 and gr_est2 have been eliminated  ### 
# 12.10.2016 (Ar): com parameter in gr.all0() is modified:   ###
#                  com=com+1  
# - 26.06.2019 (Arzideh): different percentiles:             ###
#                          pc1, pc2,PC1,PC2                  ###
################################################################
################################################################
gr.est0<-
function (data=te$dbm,est=est1,lam=est$y0[1,8],q1=est$y0[1,1],q2=est$y0[1,2],sig=est$y0[1,7],mu=est$y0[1,6],p=est$y0[1,5],
x12=0.001,x01=min(data),x02=max(data),over=40000,step0=0.1,com=1,er1=0.5*(10^(-com)),
nb=4096,c1=0.8,c2=1,sc1=1.2,labx="non-transformed values (U/l)",bw,low,high,main="estimation",lam2=mind,
pc1=0.025,pc2=0.05,PC1=0.975,PC2=0.95) 
{
# gr.est : written at 12th August 2008
# display the estimated distributions for pathological and non-pathological values and
# calculates  Enscheidungsgrenze (DL) and fase positive and false negative prob.
########################################## 
res<-gr.all0(d=data,lam=lam,q1=q1,q2=q2,sig=sig,mu=mu,p=p,lam2=lam2,
x01=x01,x02=x02,x12=x12,step0=step0,bwd=bw,over=over,er1=er1,nb=nb,
sc1=sc1,c2=c2,c1=c1,com=com+1,labx=labx,low=low,high=high,main=main,pc1=pc1,pc2=pc2,PC1=PC1,PC2=PC2)
return(list(P=res$P,mu=res$mu,sig=res$sig,T2=res$T2,DL975=res$DL975,DL95=res$DL95,myplot=res$myplot,DL25=res$DL25))
}



#cat("R function KD_tnorm.R loaded....\n")
##########################################
KD.tnorm<-function (data,mu,sig,T1,T2) 
{
# KD.tnorm evaluates the kolmogorov- distance of distr. of empirical data (cdf) from 
# estimated truncatedNormal-distribution between T1 and T2.
m<-length(data)
data1<-sort(data)
p1<-pnorm(T1,mu,sig)
p2<-pnorm(T2,mu,sig)
p0<-p2-p1
D<-c(rep(0,m));
i<-1
while (i<=m)
         {
        fn1<-i/m
        fn0<-(i-1)/m
        f0<-(pnorm(data1[i],mu,sig)-p1)/p0
        D1<-abs(fn0-f0);
        D2<-abs(f0-fn1);
        D[i]<-max(D1,D2);
        i<-i+1
        }
index<-which.max(D)
v<-data1[index]
j<-(index)/m
K<-max(D[1:m])
return(list(K=K,v=v,j=j))
}



#cat("R function KD_tnorm1.R loaded....\n")
############################################
KD.tnorm1<-function (data,p,mu,sig,T1,crit) 
{
# KD.tnorm1 evaluates the modified kolmogorov- distance of distr. of empirical data (cdf) from 
# estimated truncated Normal distribution in interval [mod,max(data)]
p1<-pnorm(T1,mu,sig)
data1<-data[data>=T1] # eliminate extreme values
data1<-sort(data1)
m<-length(data1[data1 > crit])
n<-length(data1)
n1<-length(data1[data1 <= crit])
D<-c(rep(0,m));
f00<-(p*(pnorm(data1[n1],mu,sig)-p1))/(1-p1);
i<-1
while (i<=m)
          {
          fn0<-(i-1)/n+f00;
          f0<-(p*(pnorm(data1[i+n1],mu,sig)-p1))/(1-p1);
          D1<-f0-fn0
          D[i]<-max(D1,0);
          i<-i+1
          }
D<-max(D[1:m])
return(D)
}

#cat("R function KD_tnorm2.R loaded....\n")
###########################################
KD.tnorm2<-function (data,p,mu,sig,crit) 
{
# KD.tnorm2 evaluates the modified kolmogorov- distance of distr. of empirical data (cdf) from 
# estimated truncated Normal distribution in interval [0,mod(data)]
data1<-sort(data)
n<-length(data1)
m<-length(data1[data1 <= crit])
D<-c(rep(0,m));
i<-1
while (i<=m)
          {
          fn0<-(i-1)/n;
          f0<-p*(pnorm(data1[i],mu,sig));
          D1<-f0-fn0
          D[i]<-max(D1,0);
          i<-i+1
          }
D<-max(D[1:m])
return(D)
}

#cat("R function kde.R loaded....\n")
#is loaded in program1
#############################################
kde<-
function (data,q1=0.025,x3=quantile(data,q1),nb=900,bw0=bw1)
{
# evaluation of mod of data and y_mod using density function  
t1<-density(x=data,bw=bw0,kernel="gaussian",from=x3,na.rm=TRUE,cut=3,n=nb)
k1<-t1$y
u1<-t1$x
a.max<-which.max(k1)
k<-k1[a.max]
u<-u1[a.max]
return(list(k=k,u=u))
}



#cat("R function kde1.R loaded....\n")
# is loaded in program1
##########################################
kde1<-
function (data,q1=0,x3=quantile(data,q1),nb=2048,b0=1E+2,bw0,mod)
{
# 22.04.09
# as kde function delivers the KDE of the data (k1 and u1 are the values of x and y),
# additionally it gives the values of (k1,u1) over the mod (k2,u2),
# thereby it is assumed that the pathological values are obtained only over the mod of the data.
# The values of (k2,u2) are used to estimate the curve of pathological values
########################################
######################################## 
t1<-density(data,bw=bw0,kernel="gaussian",from=x3,na.rm=TRUE,cut=3,n=nb)
k1<-t1$y
u1<-t1$x
u2<-u1[u1>=mod]
k2<-k1[u1>=mod]
m2<-length(k2)
u3<-u1[u1< mod]
k3<-k1[u1< mod]
m3<-length(k3)
bw<-t1$bw
return(list(bw1=bw,k1=k1,u1=u1,k2=k2,u2=u2,m2=m2,k3=k3,u3=u3,m3=m3))
}


# kalkuliert die RLs:
# library(msm)
############################
#cat("R function RLs.R loaded....\n")
# 12.10.2016 (Arzideh): parameter com was added. Estimated RLs are now rounded at com+1 decimals.
############################
RLs<-function(est,m=est$y0[1,6],s=est$y0[1,7],lam=est$y0[1,8],pc1=0.025,pc2=0.05,com=com)
{
# estimated RLs with minimum test-statistic
if (lam==0)
           {
             # RLs for log-Normal
             left0<-qnorm(pc1, mean=m, sd=s)
             right0<-qnorm(1-pc1, mean=m, sd=s)
             left0<-exp(left0)
             right0<-exp(right0)
             left1<-qnorm(pc2, mean=m, sd=s)
             right1<-qnorm(1-pc2, mean=m, sd=s)
             left1<-exp(left1)
             right1<-exp(right1)
            }
else
           {
            # RLs for Power-Normal
            left0<-qtnorm(pc1, mean=m, sd=s, lower=-1/lam, upper=Inf)
            right0<-qtnorm(1-pc1, mean=m, sd=s, lower=-1/lam, upper=Inf)
            left0<-(left0*lam+1)^(1/lam)
            right0<-(right0*lam+1)^(1/lam)
            left1<-qtnorm(pc2, mean=m, sd=s, lower=-1/lam, upper=Inf)
            right1<-qtnorm(1-pc2, mean=m, sd=s, lower=-1/lam, upper=Inf)
            left1<-(left1*lam+1)^(1/lam)
            right1<-(right1*lam+1)^(1/lam)
           }
return(list("RLL1"=round(left0,com+1),"RLL2"=round(left1,com+1),"RLR1"=round(right0,com+1),"RLR2"=round(right1,com+1)))
}


#cat("R function T2point.R loaded....\n")
#########################################
Tpoint<-
function (data,mod=1,com=2,c1=3,c2=3,TL1=min(data),TL2=mod-c1*(10^(-com)),TR1=mod+c2*(10^(-com)),TR2=max(data)) 
{
#  Tpoint function calculates truncation points of data
#  truncation points are determined at the left and right sides of mod of the data
#  at the left side in interval: TL=[TL1,TL2] and at the right side in interval: TR=[TR1,TR2]
#  mod, c1,c2 and com: variables to calculate TL2 and TR1, upper limit of TL and lower limit of TR respectively 
#  mod: mod of data
#  com: accuracy of data values; e.g. for data values 2.11, 1.67 is com=2 and for 144,150 is com=0
#  c1 and c2 define how far from mod TL2 and TR1 are searched
# e.g. for creatinin with data given with two decimals (0.94 or 1.03) and mod=0.88 and for c1=c2=10, TL1, TL2, TR1 and TR2 can be calculated as:
# TL1=min(data) , TR2=max(data),
# TL2= 0.88-10*(10^-2)=0.78, and TR1= 0.88+10*(10^-2)=0.98, and thereby:
# TL=[TL1,0.78], TR=[0.98,TR2]
#  lst modification: 26.09.2019: to ensure to have some data in defined interval: if(TR1 > TR2) TR2<-TR1+(10)**(-com) and  if(TL2 < TL1) TL1<-TL2-(10)**(-com)
################################################################  
################################################################ 
 n<-length(data)
 data1<-sort(data)
 h<-T<-rep(0,n)
 i<-1
 while(i<=n-1)
 {
  if (data1[i]!=data1[i+1])
   {
   h[i]<-(i)/n
   T[i]<-data1[i]
   }
  else 
   {
   h[i]<-0
   T[i]<-0
   }
  i<-i+1
 }
 h[n]<-1
 T[n]<-data1[n]
 h1<-h[h!=0]
 T<-T[T!=0]
############################
# this part ensures that some data between TR1 and TR2 exist. If TR1 and TR2 are too near (e.g equal) 
# we make TR2 larger: 
 if(TR1 > TR2) TR2<-TR1+(10)**(-com)   
 dataTR<-data1[data1>=TR1 & data1<=TR2]
 m1<-length(dataTR)
 j<-1
 while (j<=100 & m1==0)
       {
        TR2<-TR2+(10^(-com))
        dataTR<-data1[data1>=TR1 & data1<=TR2]
        m1<-length(dataTR)
        j<-j+1
        }
##########################
# this part ensures that some data between TL1 and TL2 exist. If TL1 and TL2 are too near (e.g equal) 
# we make TL2 larger:
 if(TL2 < TL1) TL1<-TL2-(10)**(-com)
 dataTL<-data1[data1>=TL1 & data1<=TL2]
 n1<-length(dataTL)
 k<-1
 while (k<=100 & n1==0)
       {
        TL2<-TL2+(10^(-com))
        dataTL<-data1[data1>=TL1 & data1<=TL2]
        n1<-length(dataTL)
        k<-k+1
        }
############################
 q2<-h1[T>=TR1 & T<=TR2]
 T2<-T[T>=TR1 & T<=TR2]
 q1<-h1[T>=TL1 & T<=TL2]
 T1<-T[T>=TL1& T<=TL2]
 if (TL1==0)
 {
 q1<-c(0,q1)
 T1<-c(0,T1)
}
 return(list(q2=q2,T2=T2,q1=q1,T1=T1))
# T1 and T2 are  vectors of truncation points and q1 and q2 the fractions of the data <=T1 and data <=T2 , respec.
}
###################################
#############################
T2point<-
function (data,mod=1,com=2,c2=3,TR1=mod+c2*(10^(-com)),TR2=max(data)) 
{
#  T2point function calculates truncation points of data
#  truncation points are determined at the right side of mod of the data in interval: TR=[TR1,TR2]
#  mod,c2 and com: variables to calculate TR1, lower limit of TR 
#  mod: mod of data
#  com: accuracy of data values; e.g. for data values 2.11, 1.67 is com=2 and for 144,150 is com=0
# c2 defines how far from mod TR1 is searched
# e.g. for creatinin with data given with two decimals (0.94 or 1.03) and mod=0.88 and for c2=10, TR1 and TR2 can be calculated as:
# , TR2=max(data),
#  TR1= 0.88+10*(10^-2)=0.98, and thereby:
#  TR=[0.98,TR2]
 n<-length(data)
 data1<-sort(data)
 h<-T<-rep(0,n)
 i<-1
 while(i<=n-1)
 {
  if (data1[i]!=data1[i+1])
   {
   h[i]<-i/n
   T[i]<-data1[i]
   }
  else 
   {
   h[i]<-0
   T[i]<-0
   }
  i<-i+1
 }
 h[n]<-1
 T[n]<-data1[n]
 h1<-h[h!=0]
 T<-T[T!=0]
#####################################
# this part ensures that some data between TR1 and TR2 exist. If TR1 and TR2 are too near (e.g equal) 
# we make TR2 larger:
 if(TR1 > TR2) TR2<-TR1+(10)**(-com)  
 dataT<-data1[data1>=TR1 & data1<=TR2]
 m1<-length(dataT)
 j<-1
 while (j<=100 & m1==0)
       {
        TR2<-TR2+(10^(-com))
        dataT<-data1[data1>=TR1 & data1<=TR2]
        m1<-length(dataT)
        j<-j+1
        } 
####################################
 q2<-h1[T>=TR1 & T<=TR2]
 T2<-T[T>=TR1 & T<=TR2] 
 return(list(q2=q2,T2=T2))
# T2 is the vector of truncation points and q2 the fraction of the data <=T2, respec.
}
#######################
T1point<-
function (data,mod=1,com=2,c2=3,TL2=mod-c2*(10^(-com)),TL1=min(data)) 
{
#  T1point function calculates truncation points of data
#  truncation points are determined at the left side of mod of the data in interval: TL=[TL1,TL2]
#  mod,c2 and com: variables to calculate TR1, lower limit of TR 
#  mod: mod of data
#  com: accuracy of data values; e.g. for data values 2.11, 1.67 is com=2 and for 144,150 is com=0
#  c2 defines how far from mod TL2 is searched
#  e.g. for creatinin with data given with two decimals (0.94 or 1.03) and mod=0.88 and for c2=10, TL1 and TL2 
#  can be calculated as:
#  TL1=min(data),
#  TL2= 0.88-10*(10^-2)=0.78, and thereby:
#  TL=[TL1,0.78]
 n<-length(data)
 data1<-sort(data)
 h<-T<-rep(0,n)
 i<-2
 while(i<=n)
 {
  if (data1[i]!=data1[i-1])
   {
   h[i]<-(i-1)/n
   T[i]<-data1[i]
   }
  else 
   {
   h[i]<-0
   T[i]<-0
   }
  i<-i+1
 }
 h[1]<-0
 T[1]<-data1[1]
 h1<-h[T!=0]
 T<-T[T!=0]
##########################
# this part ensures that some data between TL1 and TL2 exist. If TL1 and TL2 are too near (e.g equal) 
# we make TL2 larger:
 if(TL2 < TL1) TL1<-TL2-(10)**(-com)
 dataTL<-data1[data1>=TL1 & data1<=TL2]
 n1<-length(dataTL)
 k<-1
 while (k<=100 & n1==0)
       {
        TL2<-TL2+(10^(-com))
        dataTL<-data1[data1>=TL1 & data1<=TL2]
        n1<-length(dataTL)
        k<-k+1
        }
########################
 q1<-h1[T>=TL1 & T<=TL2]
 T1<-T[T>=TL1 & T<=TL2]
 return(list(q1=q1,T1=T1))
# T1 is the vector of truncation points and q1 the fraction of the data > = T1, respec.
}


#cat("R function TNorm.R loaded....\n")
# last change: 24.11.2021: any(class(iHess) == 'try-error' [see below]!
#########################################
TNorm<-
function(d,q1=0.0,T1=quantile(d,q1,na.rm =T),q2,T2=quantile(d,q2,na.rm =T),I=100,epsilon1=1E-4,
mu,sig,er1=0.0)
 {
  # Estimation of parameters of a truncated normal distribution
  # q1 and q2 are left and right quantiles, at which the data is truncated.
  # mu and sig are initial values for estimation.
  q1<-(length(d[d<=T1]))/length(d)
  q2<-(length(d[d<=T2]))/length(d)
  d.trunc<-d[d>T1 & d<T2]
  n<-length(d.trunc); m<-mean(d.trunc); s<-sd(d.trunc); i<-1; mu0<-mu; sig0<-sig; l.old<-1; 
  psi1<-(T1-mu)/sig; psi2<-(T2-mu)/sig;
  Phi1<-pnorm(psi1); Phi2<-pnorm(psi2); Phi<-Phi2-Phi1;
  theta.new<-t(t(c(mu,sig)));
  Diff<-t(t(c(0,0)));
  Ide<-diag(rep(1,2))
  differential<-1
  Hess<-iHess<-matrix(rep(0,4),ncol=2); 
  l.new<-(-n)*(log(Phi))-(n)*(log(sig))-(n/2)*(log(2*pi))-(1/(2*(sig^2)))*(sum((d.trunc-mu)^2))
  cr<-1
  while( i<=I && (differential > epsilon1))
          {#while
          phi1<-dnorm(psi1); phi2<-dnorm(psi2);
          Q1<-phi1/Phi; Q2<-phi2/Phi;
          P1<-Q1*(Q1-psi1);   P2<-(-Q2)*(Q2+psi2);
          la1<-psi1*P1+Q1; la2<-psi2*P2+Q2;
          et1<-psi1*(la1+Q1); et2<-psi2*(la2+Q2);
          jadid<-(-2*psi1*psi2*Q1*Q2)
          n1<-n/sig^2;
         Diff1<-n1*(m-mu-sig*(Q1-Q2));
         Diff2<-n1*((s^2+(m-mu)^2)/sig-sig*(1+psi1*Q1-psi2*Q2));
         Diff<-t(t(c(Diff1,Diff2)));
         Hess[1,1]<-(-n1)*(1-P1+P2);
         Hess[1,2]<-Hess[2,1]<-(-n1)*(2*(m-mu)/sig-la1+la2);
         Hess[2,2]<-(-n1)*(3*(s^2+(m-mu)^2)/sig^2-1-et1+et2-jadid);
         iHess<-try(qr.solve(Hess),TRUE);
         if ((any(is.na(Hess))) | abs(det(Hess)) < 1E-4 | any(class(iHess) == 'try-error'))
            {#if1
                                   l.old<-1;
                                   sig<-mu<-Phi<-l.new<-frac<-differential<-cr<-1E+5
                                   iHess<-matrix(rep(1,4),ncol=2) 
                            i<-I+1
         } else {  
         theta.old<-theta.new;
         theta.new<-theta.old-(iHess)%*%(Diff)  
         mu<-theta.new[1,1];
         sig<-theta.new[2,1];
         mu.old<-theta.old[1,1];
         sig.old<-theta.old[2,1];
         di1<-abs(mu-mu.old)
         di2<-abs(sig-sig.old)
         cr<-max(di1,di2)
         differential<-sqrt(Diff1^2+Diff2^2)
         l.old<-l.new; 
         psi1<-(T1-mu)/sig; psi2<-(T2-mu)/sig;
        Phi1<-pnorm(psi1); Phi2<-pnorm(psi2);
        Phi<-Phi2-Phi1;
        #if (sig > 1E-20 & Phi > 1E-20 & mu > T1)
        if (sig > 1E-20 & Phi > 1E-20 )
        {#if2                            
        l.new<-(-n)*(log(Phi))-(n)*(log(sig))-(n/2)*(log(2*pi))-(1/(2*(sig^2)))*(sum((d.trunc-mu)^2));
        frac<-min((q2-q1)/Phi,1);
        i<-i+1
         }#if2
        else 
                           {#else2
                           l.old<-1;
                           sig<-mu<-Phi<-l.new<-frac<-differential<-cr<-1E+5
                            i<-I+1
                            }#else2
       }#else1
       }#while
 return(list(q1=q1,q2=q2,T1=T1,T2=T2,I=i-1,mu=mu,sig=sig,Phi=Phi,frac=frac,l=l.new,
Diff=Diff,cov=-iHess,differential=differential))
}




tml <- function(d, pathright) {
  
  ### variables from parameters.txt
  Q1 <- 0.10
  Q2 <- 0.90
  dec.num <- 1
  minsize <- 40
  dif1 <- 1
  nen <- 1
  quant <- 3
  pathol <- "AUTO"
  q10 <- quantile(d, probs = 0.1, na.rm=TRUE)
  q90 <- quantile(d, probs = 0.9, na.rm=TRUE)
  x1 <- q10 - (q90 - q10) / 1.3
  x2 <- q90 + (q90 - q10) / 1.3
  ci <- FALSE
  nrep <- 500
  model <- "PN"
  alpha1 <- 5
  l0 <- 0
  n <- 10
  lstep <- 0.1
  nstep <- 21
  minsize <- 4000
  ci_level <- 0.9
  low_RL <- 0.025
  upp_RL <- 0.975
  meth.all <- 0
  
  ### KDE
  if (length(d) > 40000) {
    bw <- bw.nrd(d)
  } else{
    bw <- try(bw.SJ(d, method = "dpi"), silent = T)
    if (class(bw) == "try-error") {
      bw <- bw.nrd(d)
    }
  }
  step0 <- (10 ^ (-dec.num))
  bw <- max(bw, (step0) / 2)
  step <- max(bw, step0)
  mod <- kde(
    data = d,
    nb = 900,
    q1 = 0,
    bw0 = bw
  )$u
  
  ### non path

if (pathright) estimate <-
    EST.PN(
      data = d,
      Q1 = Q1,
      Q2 = Q2,
      dec = dec.num,
      mod = mod,
      l0 = l0,
      n = n,
      dif1 = dif1,
      nen = nen,
      lstep = lstep,
      nstep = nstep
    )
else estimate <-
    EST.PN.L(
      data = d,
      Q1 = Q1,
      Q2 = Q2,
      dec = dec.num,
      mod = mod,
      l0 = l0,
      n = n,
      dif1 = dif1,
      nen = nen,
      lstep = lstep,
      nstep = nstep
    )

  esta <- estimate$est
  ta <- RLs(
    est = esta,
    com = dec.num,
    pc1 = low_RL,
    pc2 = upp_RL
  )
  TLa <- esta$y0[1, 3]
  TR1a <- esta$y0[1, 4]
  
  ### go
  return(
    gr.est0(
      data = d,
      est = esta,
      x02 = x2,
      x01 = x1,
      step0 = step,
      bw = bw,
      com = dec.num,
      labx = NULL,
      low = TLa,
      high = TR1a,
      main = NULL,
      lam2 = 0,
      pc1 = low_RL,
      PC1 = upp_RL
    ) # $myplot
  )
}
