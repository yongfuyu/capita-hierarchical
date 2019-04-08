model_string<-"
model{
for(i in 1:N){
N_IPD[i]~dpois(lambda[i])
log(lambda[i])<- (beta[1,st.index[i]] +beta[2,st.index[i]]*vax[i] +log(N_POP[i]))
}

for(j in 1:2){
mu[j]~dnorm(0,1e-4)
tau[j]<-1/sd1[j]^2
sd1[j]~dunif(0,40)
for(k in 1:13){
beta[j,k]~dnorm(mu[j], tau[j])
}
}

}
"

##############################################################
#Model Fitting
##############################################################
model_jags<-jags.model(textConnection(model_string),n.chains=2,
                       data=list('N' = nrow(d1),
                                 'st.index'=d1$st.index,
                                 'N_IPD'=d1$n_s11_pp,
                                 'N_POP'=d1$pop,
                                 'vax'=d1$vax
                       )) 

update(model_jags, 
       n.iter=10000) 

posterior_samples<-coda.samples(model_jags, 
                                variable.names=c("mu",
                                                 "sd1",
                                                 "beta"),
                                thin=1,
                                n.iter=10000)

post1.summary<-summary(posterior_samples)
plot(posterior_samples, 
     ask=TRUE)
coefs<-post1.summary[['quantiles']][,c('2.5%','50%','97.5%')] 
st.slp.rows<-grep('beta[2,', row.names(coefs), fixed=T)
st.labs<-as.character(unique(d1$st))
post.all<-do.call(rbind,posterior_samples)
plot(post.all[,'beta[2,3]'], type='l') 
hist(post.all[,'beta[2,3]']) 
prob.effect.st.u0<-apply(post.all[,st.slp.rows],2,function(x) mean(x<0))
prob.effect.st.u0
log.rr.st<-t(apply(post.all[,st.slp.rows],2,quantile, probs=c(0.025,0.5,0.975)))

global.slp.rows<-grep('mu[2', row.names(coefs), fixed=T)
overall.VE<-100*round( 1 - exp(coefs[global.slp.rows,]) ,2) 
#st.VE<- cbind(as.character(unique(d1$st)),round( 1 - exp(coefs[st.slp.rows,]) ,2))
st.VE<- 100*round( 1 - exp(coefs[st.slp.rows,]) ,2)

yrange<-range(st.VE)

combined.ve<-rbind(overall.VE, st.VE)

#install.packages('rmeta')
library(rmeta)
summary_data <- 
  structure(list(
    mean  = c(NA, NA,st.VE[,'50%'],NA,overall.VE['50%']), 
    lower = c(NA, NA,st.VE[,'97.5%'],NA,overall.VE['97.5%']),
    upper = c(NA, NA,st.VE[,'2.5%'],NA,overall.VE['2.5%']),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -11L), 
    class = "data.frame"))

tabletext<-cbind(
  c("", "Serotype", st.labs, NA, "All"),
  c("Vaccine", "(N=42240)", d1$n_s11_pp[d1$vax==1], NA, sum(d1$n_s11_pp[d1$vax==1])),
  c("Control", "(N=42256)",  d1$n_s11_pp[d1$vax==0], NA, sum(d1$n_s11_pp[d1$vax==0])),
  c("", "VE", paste0(st.VE[,'50%'],'%'), NA, paste0(overall.VE['50%'],'%') ))

rmeta::forestplot(tabletext, 
                  mean=summary_data$mean,lower=summary_data$lower,
                  upper=summary_data$upper,
                  new_page = F,
                  #is.summary=c(rep(F, 16), T),
                  is.summary=F,
                  #clip=c(-100,500), 
                  xlog=F, 
                  boxsize=0.5)
