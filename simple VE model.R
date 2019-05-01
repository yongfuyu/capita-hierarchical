simp.mod.func<-function(){
model_string<-"
model{
for(i in 1:N){
N_IPD[i]~dpois(lambda[i])
log(lambda[i])<- (mu[1] +mu[2]*vax[i] +log(N_POP[i]))
}

for(j in 1:2){
mu[j]~dnorm(0,1e-4)
tau[j]<-1/sd1[j]^2
sd1[j]~dunif(0,40)
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
                                 'N_IPD'=d1$n_s11_pp,
                                 'N_POP'=d1$pop,
                                 'vax'=d1$vax
                       )) 

update(model_jags, 
       n.iter=10000) 

posterior_samples<-coda.samples(model_jags, 
                                variable.names=c("mu"),
                                thin=1,
                                n.iter=50000)

post1.summary<-summary(posterior_samples)
log.rr<-post1.summary[[2]]['mu[2]',c('2.5%', '50%', '97.5%')]
ve.simple<- round(100*(1-exp(log.rr)), 2) 
return(ve.simple)
}