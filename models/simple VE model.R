simp.mod.func<-function(){
model_string<-"
model{
for(i in 1:N){
N_IPD[i]~dpois(lambda[i])
log(lambda[i])<- (beta[1,st.index[i]] +mu*vax[i] +log(N_POP[i]))
}

for(k in 1:13){
beta[1,k]~dnorm(0,1e-4)
}

mu~dnorm(0,1e-4)

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
                                 'N_IPD'=d1[,outcome.var],
                                 'st.index'=d1$st.index,
                                 'N_POP'=d1$pop,
                                 'vax'=d1$vax
                       )) 

update(model_jags, 
       n.iter=10000) 

posterior_samples<-coda.samples(model_jags, 
                                variable.names=c("mu"),
                                thin=1,
                                n.iter=50000)
ci<-p.interval(posterior_samples[[1]],
                   HPD=TRUE,
                   MM=TRUE) 
#post1.summary<-summary(posterior_samples)
#log.rr<-post1.summary[[2]][c('2.5%', '50%', '97.5%')]
post_means<-colMeans(posterior_samples[[1]])

ve.simple<- round(100*(1-exp(c(post_means, ci))), 2) 
return(ve.simple)
}