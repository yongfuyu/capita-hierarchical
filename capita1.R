###REPORTED EFFICACY
#Vaccine-type CAP PP: 45.6(21.8, 62.5) (Tabl S10): n_s10_pp
#Vaccine-type CAP ITT: 37.7(14.3, 55.1)
#Vaccine-type Non-bacteremic CAP PP: 45.0(14.2, 65.3) Table S11: n_s11_pp
#Vaccine-type Non-bacteremic CAP ITT: 41.1(12.7, 60.7)
library(rjags)
library(rmeta)

d1<-read.csv('./capita.csv')
d1$st.index<-rep(1:13, each=2)
st.n<- d1$n_s11_pp[d1$vax==0] +d1$n_s11_pp[d1$vax==1]
mod1<-glm(n_s11_pp~ st +vax+ offset(log(pop)), data=d1, family='poisson')
vax.coef<-coef(summary(mod1))['vax','Estimate']
vax.coef.se<-coef(summary(mod1))['vax','Std. Error']
VE<- 1-round(exp(vax.coef),2)
VE.lcl<-round(1-exp(vax.coef-1.96*vax.coef.se),2)
VE.ucl<-round(1-exp(vax.coef+1.96*vax.coef.se),2)
paste0(VE, '(',VE.lcl,', ', VE.ucl,')')

#####################################################################################
#Simple VE model without serotype
#####################################################################################
source('./models/simple VE model.R')
ve.simple<-simp.mod.func()
ve.simple

#####################################################################################
#Non-Hierarchical Model
#####################################################################################
source('./models/non-hierarchical model.R')
nonhier1<-non_hierarchical_func()
nonhier.st.eff<-nonhier1$st.VE
pdf('./results/non_hierarchical.pdf', width=14, height=7)
rmeta::forestplot(nonhier1$tabletext, 
                  mean=nonhier1$summary_data$mean,
                  lower=signif(nonhier1$summary_data$lower,digits=2),
                  upper=signif(nonhier1$summary_data$upper, digits=2),
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
hier.st.eff<-hier1$st.VE

pdf('./results/hierarchical.pdf', width=8, height=7)
rmeta::forestplot(hier1$tabletext, 
                         mean=hier1$summary_data$mean,lower=hier1$summary_data$lower,
                         upper=hier1$summary_data$upper,
                         new_page = F,
                         #is.summary=c(rep(F, 16), T),
                         is.summary=F,
                         clip=c(-50,100), 
                         xlog=F, 
                        align='l',
                         boxsize=0.5)
dev.off()



##############################
#HIERARCHICAL MIXTURE MODEL
##############################
source('./models/hierarchical mixture model.R')
hier.mix1<-hierarchical.mix.mod.func()
hier.mix.st.eff<-hier.mix1$st.VE

pdf('./results/hierarchical mixture.pdf', width=8, height=7)
rmeta::forestplot(hier1$tabletext, 
                  mean=hier.mix1$summary_data$mean,lower=hier.mix1$summary_data$lower,
                  upper=hier.mix1$summary_data$upper,
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
pdf('./results/hierarchical vs non hierarchical xy.pdf', width=8, height=8)
par(mfrow=c(1,1))
plot( nonhier.st.eff[,2],hier.st.eff[,2], ylim=c(-100,105), 
     xlim=c(-100,105), bty='l', ylab='Hierarchical', xlab='non-hierarchical')
#arrows(x0=hier.st.eff[,2], x1=hier.st.eff[,2], y0=nonhier.st.eff[,1],y1=nonhier.st.eff[,3], length=0 )
#arrows(x0=hier.st.eff[,1], x1=hier.st.eff[,3], y0=nonhier.st.eff[,2],y1=nonhier.st.eff[,2], length=0 )
text( nonhier.st.eff[,2],hier.st.eff[,2], unique(d1$st), adj=c(1,1), cex=0.5)
abline(h=0,v=0, col='gray', lty=2)
#abline(h=60, v=50)
dev.off()

pdf('./results/hierarchical vs non hierarchical shrinkage.pdf', width=8, height=6)
plot.nonhier.st.eff <- nonhier.st.eff[,2]
plot.nonhier.st.eff[plot.nonhier.st.eff< -10000]<-NA
plot( plot.nonhier.st.eff,rep(2,13), ylim=c(1,2.1),yaxt='n',ylab='',xlab='Efficacy', xlim=c(-200,100),col='white', bty='n')
symbols(plot.nonhier.st.eff,rep(2,13),st.n, cex=0.2, inches=0.1, add=T)
#points(hier.st.eff[,2],rep(1,13))
symbols(hier.st.eff[,2],rep(1,13),st.n, cex=0.2, inches=0.1, add=T)
arrows(x0=plot.nonhier.st.eff , x1=hier.st.eff[,2], y0=2, y1=1, length=0, lty=2, col='gray')
text(-145,2, 'Non-hierarchical', pos=2)
text(-145,1, 'Hierarchical', pos=2)
abline(v=0) 
dev.off()

