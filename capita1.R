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
#Non-Hierarchical Model
#####################################################################################
source('./non-hierarchical model.R')

#####################################################################################
#Hierarchical Model
#####################################################################################
source('./hierarchical model.R')
