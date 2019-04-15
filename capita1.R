###REPORTED EFFICACY
#Vaccine-type CAP PP: 45.6(21.8, 62.5) (Tabl S10): n_s10_pp
#Vaccine-type CAP ITT: 37.7(14.3, 55.1)
#Vaccine-type Non-bacteremic CAP PP: 45.0(14.2, 65.3) Table S11: n_s11_pp
#Vaccine-type Non-bacteremic CAP ITT: 41.1(12.7, 60.7)

library(rjags)
library(rmeta)

d1<-read.csv('./capita.csv')
d1$st.index<-rep(1:13, each=2)
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
source('./simple VE model.R')
ve.simple<-simp.mod.func()
ve.simple

#####################################################################################
#Non-Hierarchical Model
#####################################################################################
source('./non-hierarchical model.R')
nonhier1<-non_hierarchical_func()
nonhier.st.eff<-nonhier1$st.VE
pdf('non_hierarchical.pdf', width=14, height=7)
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
source('./hierarchical model.R')
hier1<-hierarchical.mod.func()
hier.st.eff<-hier1$st.VE
pdf('hierarchical.pdf', width=8, height=7)

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
par(mfrow=c(1,1))
plot(hier.st.eff[,2], nonhier.st.eff[,2], ylim=c(-100,105), 
     xlim=c(-100,105), bty='l', ylab='Non-hierarchical', xlab='Hierarchical')
#arrows(x0=hier.st.eff[,2], x1=hier.st.eff[,2], y0=nonhier.st.eff[,1],y1=nonhier.st.eff[,3], length=0 )
#arrows(x0=hier.st.eff[,1], x1=hier.st.eff[,3], y0=nonhier.st.eff[,2],y1=nonhier.st.eff[,2], length=0 )
abline(h=0,v=0, col='gray', lty=2)
