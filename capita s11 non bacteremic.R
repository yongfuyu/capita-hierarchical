###REPORTED EFFICACY
#Vaccine-type CAP PP: 45.6(21.8, 62.5) (Tabl S10): n_s10_pp
#Vaccine-type CAP ITT: 37.7(14.3, 55.1)
#Vaccine-type Non-bacteremic CAP PP: 45.0(14.2, 65.3) Table S11: n_s11_pp
#Vaccine-type Non-bacteremic CAP ITT: 41.1(12.7, 60.7)
library(rjags)
library(rmeta)
library(HDInterval)

#Name of the outcome variable
########################
outcome.var<- 'n_s11_pp'
########################

d1<-read.csv('./capita.csv')
d1$st.index<-rep(1:13, each=2)
st.n<- d1[,outcome.var][d1$vax==0] +d1[,outcome.var][d1$vax==1]

#######################################################
#Preparing Data
#######################################################
m<-c(42256, 42240)

vax<-c(0, 1)

y<-matrix(0, nrow=2, ncol=14)
y[1, 1:13]<-d1[d1$vax==0,outcome.var]
y[2, 1:13]<-d1[d1$vax==1,outcome.var]

y[1,14]<-m[1] - sum(y[1, 1:13])
y[2,14]<-m[2] - sum(y[2, 1:13])

#####################################################################################
# Model with same VE for all serotypes
#####################################################################################
source('./models/simple VE model.R')
ve.non.varying<-simp.mod.func()
ve.non.varying

#####################################################################################
#Non-Hierarchical Model
#####################################################################################
source('./models/non-hierarchical model.R')
nonhier1<-non_hierarchical_func()
nonhier.st.eff<-as.numeric(nonhier1$st.VE[,1])
tiff('./results/non_hierarchical.tif', width=8, height=7,res=200,units = 'in')
rmeta::forestplot(nonhier1$tabletext, 
                  mean=as.numeric(nonhier1$summary_data$mean),
                  lower=as.numeric(nonhier1$summary_data$lower),
                  upper=as.numeric(nonhier1$summary_data$upper),
                  new_page = T,
                  clip=c(-50,100),
                  #is.summary=c(rep(F, 16), T),
                  is.summary=F,
                  align='l',
                  xlog=F, 
                  boxsize=0.5)
dev.off()
#####################################################################################
#Hierarchical Model
#####################################################################################
source('./models/hierarchical model.R')
hier1<-hierarchical.mod.func()
hier.st.eff<-as.numeric(hier1$st.VE[,1])

tiff('./results/hierarchical.tif', width=8, height=7,res=200,units = 'in')
rmeta::forestplot(hier1$tabletext, 
                         mean=as.numeric(hier1$summary_data$mean),
                          lower=as.numeric(hier1$summary_data$lower),
                         upper=as.numeric(hier1$summary_data$upper),
                         new_page = F, 
                         #is.summary=c(rep(F, 16), T),
                         is.summary=F,
                         clip=c(-50,100), 
                         xlog=F, 
                        align='l',
                         boxsize=0.5)
dev.off()



###############################
##COMPARE HIERARCHICAL VS NON-HIERARCHICAL
tiff('./results/hierarchical vs non hierarchical xy.tif', width=8, height=8,res=200,units = 'in')
par(mfrow=c(1,1))
plot( nonhier.st.eff,hier.st.eff, ylim=c(-100,105), 
      xlim=c(-100,105), bty='l', ylab='Hierarchical', xlab='non-hierarchical')
#arrows(x0=hier.st.eff[,2], x1=hier.st.eff[,2], y0=nonhier.st.eff[,1],y1=nonhier.st.eff[,3], length=0 )
#arrows(x0=hier.st.eff[,1], x1=hier.st.eff[,3], y0=nonhier.st.eff[,2],y1=nonhier.st.eff[,2], length=0 )
text( nonhier.st.eff,hier.st.eff, unique(d1$st), adj=c(1,1), cex=0.5)
abline(h=0,v=0, col='gray', lty=2)
#abline(h=60, v=50)
dev.off()

tiff('./results/hierarchical vs non hierarchical shrinkage.tif', width=8, height=6,res=200,units = 'in')
plot.nonhier.st.eff <- nonhier.st.eff
plot.nonhier.st.eff[plot.nonhier.st.eff< -10000]<-NA
plot( plot.nonhier.st.eff,rep(2,13), ylim=c(1,2.1),yaxt='n',ylab='',xlab='Vaccine Efficacy', xlim=c(-200,100),col='white', bty='n')
symbols(plot.nonhier.st.eff,rep(2,13),sqrt(st.n/pi), cex=0.2, inches=0.1, add=T)
#points(hier.st.eff[,2],rep(1,13))
symbols(hier.st.eff,rep(1,13),sqrt(st.n/pi), cex=0.2, inches=0.1, add=T)
arrows(x0=plot.nonhier.st.eff , x1=hier.st.eff, y0=2, y1=1, length=0, lty=2, col='gray')
text(-145,2, 'Non-hierarchical', pos=2)
jitter.plot.nonhier.st.eff<-plot.nonhier.st.eff
jitter.y<-rep(2.01, length(jitter.plot.nonhier.st.eff))
jitter.y[unique(d1$st)=='23F']<-2.04
jitter.y[unique(d1$st)=='19F']<-1.98
jitter.y[unique(d1$st)=='5']<-2.04
jitter.y[unique(d1$st)=='19A']<-2.06
jitter.y[unique(d1$st)=='3']<-2.03
text(jitter.plot.nonhier.st.eff,jitter.y,  unique(d1$st), pos=3, xpd=NA, cex=0.75)
text(-145,1, 'Hierarchical', pos=2)
abline(v=0) 
dev.off()



#######################################################################################################################
#SENSITIVITY ANALYSES: ALTERNATIVE PRIOR STRUCTURES
#######################################################################################################################
#####################################################################################
#Hierarchical Model with inverse gamma priors on the shrinkage parameters
#####################################################################################
source('./models/hierarchical model-IG.R')
hier2<-hierarchical.ig.mod.func()
hier.st.eff<-hier2$st.VE

tiff('./results/hierarchical-ig.tif', width=8, height=7,res=200,units = 'in')
rmeta::forestplot(hier2$tabletext, 
                  mean=as.numeric(hier2$summary_data$mean),
                  lower=as.numeric(hier2$summary_data$lower),
                  upper=as.numeric(hier2$summary_data$upper),
                  new_page = F,
                  #is.summary=c(rep(F, 16), T),
                  is.summary=F,
                  clip=c(-50,100), 
                  xlog=F, 
                  align='l',
                  boxsize=0.5)
dev.off()

###########################################################################################
###DIRICHLET PRIOR, NON_HIERARCHICAL 
source('./models/non-hierarchical model-Dirichlet.R')
nonhier2<-non_hierarchical_dir_func()
nonhier.st.eff2<-nonhier2$st.VE
tiff('./results/non_hierarchical-dirichlet.tif', width=8, height=7,res=200,units = 'in')
rmeta::forestplot(nonhier2$tabletext, 
                  mean=as.numeric(nonhier2$summary_data$mean),
                  lower=as.numeric(nonhier2$summary_data$lower),
                  upper=as.numeric(nonhier2$summary_data$upper),
                  new_page = T,
                  clip=c(-50,100),
                  #is.summary=c(rep(F, 16), T),
                  is.summary=F,
                  align='l',
                  xlog=F, 
                  boxsize=0.5)
dev.off()

