non_hierarchical_func<-function(){
model_string<-"
model{
for(i in 1:N){
N_IPD[i]~dpois(lambda[i])
log(lambda[i])<- (beta[1,st.index[i]] +beta[2,st.index[i]]*vax[i] +log(N_POP[i]))
}

for(j in 1:2){
for(k in 1:13){
beta[j,k]~dnorm(0, 1e-4)
}
}

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
                                variable.names=c("beta"),
                                thin=1,
                                n.iter=50000)

post1.summary<-summary(posterior_samples)
#plot(posterior_samples, 
#     ask=TRUE)
coefs<-post1.summary[['quantiles']][,c('2.5%','50%','97.5%')] 
st.slp.rows<-grep('beta[2,', row.names(coefs), fixed=T)
st.labs<-as.character(unique(d1$st))
post.all<-do.call(rbind,posterior_samples)
#plot(post.all[,'beta[2,3]'], type='l') 
#hist(post.all[,'beta[2,3]']) 
prob.effect.st.u0<-apply(post.all[,st.slp.rows],2,function(x) mean(x<0))
prob.effect.st.u0
log.rr.st<-t(apply(post.all[,st.slp.rows],2,quantile, probs=c(0.025,0.5,0.975)))

st.VE<- 100*round( 1 - exp(coefs[st.slp.rows,]) ,2)

yrange<-range(st.VE)

#install.packages('rmeta')
library(rmeta)
summary_data <- 
  structure(list(
    mean  = c(NA, NA,st.VE[,'50%']), 
    lower = c(NA, NA,st.VE[,'97.5%']),
    upper = c(NA, NA,st.VE[,'2.5%']),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -11L), 
    class = "data.frame"))

if (outcome.var=='n_s11_pp'){
  st.VE[8,]<-NA
}
tabletext<-cbind(
  c("", "Serotype", st.labs),
  c("Vaccine", "(N=42240)", d1[,outcome.var][d1$vax==1]),
  c("Control", "(N=42256)",  d1[,outcome.var][d1$vax==0]),
  c("", "VE", paste0(st.VE[,'50%'],'% (', st.VE[,'97.5%'],'%, ',st.VE[,'2.5%'] ,'%' ,')'))
  ) 
    
  
res.list<-list('tabletext'=tabletext, 'summary_data'=summary_data,'st.VE'=st.VE)
return(res.list)

}