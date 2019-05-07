hierarchical.mix.mod.func<-function(){
model_string<-"
model{
for(i in 1:N){
N_IPD[i]~dpois(lambda[i])
log(lambda[i])<- (beta[1,st.index[i]] +beta[2,st.index[i]]*vax[i] +log(N_POP[i]))
}

#Intercept
  mu.int[1] ~ dnorm(0, 0.0001)
  tau.int<-1/sd1.int^2
  sd1.int~dunif(0,40)
  for(k in 1:13){
    beta[1,k]~dnorm(mu.int, tau.int)
  }

#Mixture mean for slope
  mu[1] ~ dnorm(0, 0.0001)
  mu[2] <- mu[1] + delta
  delta ~ dunif(0, 100)
  tau<- 1/(sd1^2)
  sd1 ~ dunif(0, 100)
  for(k in 1:13){
    beta[2,k] ~ dnorm(mu[(g[k] + 1)], tau)
    g[k] ~ dbin(p, 1)      
  }
  

p ~ dunif(0, 1)

}
"

##############################################################
#Model Fitting
##############################################################
jags.inits1 <- function(){
  list(".RNG.seed"=c(1), ".RNG.name"='base::Wichmann-Hill')
}
model_jags<-jags.model(textConnection(model_string),n.chains=2,
                       inits = jags.inits1,
                       data=list('N' = nrow(d1),
                                 'st.index'=d1$st.index,
                                 'N_IPD'=d1[,outcome.var],
                                 'N_POP'=d1$pop,
                                 'vax'=d1$vax
                       )) 

update(model_jags, 
       n.iter=10000) 

posterior_samples<-coda.samples(model_jags, 
                                variable.names=c("mu",
                                                 "sd1",
                                                 "beta", 'delta','g', 'mu.int'),
                                thin=1,
                                n.iter=50000)

post1.summary<-summary(posterior_samples)
#plot(posterior_samples,      ask=TRUE)
coefs<-post1.summary[['quantiles']][,c('2.5%','50%','97.5%')] 
st.slp.rows<-grep('beta[2,', row.names(coefs), fixed=T)
st.labs<-as.character(unique(d1$st))
post.all<-do.call(rbind,posterior_samples)
#plot(post.all[,'beta[2,3]'], type='l') 
#hist(post.all[,'beta[2,3]']) 
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
  c("Vaccine", "(N=42240)", d1[,outcome.var][d1$vax==1], NA, sum(d1[,outcome.var][d1$vax==1])),
  c("Control", "(N=42256)",  d1[,outcome.var][d1$vax==0], NA, sum(d1[,outcome.var][d1$vax==0])),
  c("", "VE", paste0(st.VE[,'50%'],'% (', st.VE[,'97.5%'],'%, ',st.VE[,'2.5%'] ,'%' ,')'), 
    NA, paste0(overall.VE['50%'],'% (', overall.VE['97.5%'],'%, ', overall.VE['2.5%'],'%',')' ) )
  
  )

res.list<-list('tabletext'=tabletext, 'summary_data'=summary_data,'overall.VE'=overall.VE,'st.VE'=st.VE)
return(res.list)

}