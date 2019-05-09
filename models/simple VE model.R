simp.mod.func<-function(){
  model_string<-"
model{
  
  for(i in 1:2){
  
  #Likelihood
  y[i, 1:14] ~ dmulti(p[i, 1:14], m[i])
  
  #Multionmial Logistic Regression    
  for(j in 1:13){
  p_temp[i,j] <- exp(beta0[j] + beta1*vax[i])
  }
  for(j in 1:13){
  p[i,j] <- p_temp[i,j]/(1 + sum(p_temp[i, 1:13]))
  }
  p[i,14] <- 1/(1 + sum(p_temp[i, 1:13]))
  
  }
  
  #Serotype-Specific Vaccine Effects
  #100*(p_{no_vax} - p_{vax})/p_{no_vax}
  for(j in 1:13){
  sero_vax_effect[j] <- 100*(1 - p[2,j]/p[1,j])
  }
  
  #Priors
  beta1 ~ dnorm(0, 0.0001)
  for(j in 1:13){
  beta0[j] ~ dnorm(0, 0.0001)
  }
  
  #Overall Vaccine Effect
  overall_vax_effect <- 100*(1 - (1 - p[2,14])/(1 - p[1,14]))
  
}
"

##############################################################
#Model Fitting
##############################################################
#jags.inits1 = parallel.seeds("base::BaseRNG", 3) # a list of lists

inits1=list(".RNG.seed"=c(123), ".RNG.name"='base::Wichmann-Hill')
inits2=list(".RNG.seed"=c(456), ".RNG.name"='base::Wichmann-Hill')
inits3=list(".RNG.seed"=c(789), ".RNG.name"='base::Wichmann-Hill')

##############################################
#Model Organization
##############################################
model_spec<-textConnection(model_string)
model_jags<-jags.model(model_spec, 
                       inits=list(inits1,inits2, inits3),
                       data=list('y' = y,
                                 'm' = m,
                                 'vax' = vax),
                       n.adapt=10000, 
                       n.chains=3)

params<-c('overall_vax_effect')
##############################################
#Posterior Sampling
##############################################
posterior_samples<-coda.samples(model_jags, 
                                params, 
                                n.iter=10000)
posterior_samples.all<-do.call(rbind,posterior_samples)

############################################################################################################################
#Posterior Inference
############################################################################################################################
post_means<-colMeans(posterior_samples.all)
  ci<-hdi(posterior_samples.all,
              credMass = 0.95)



#post1.summary<-summary(posterior_samples)
#log.rr<-post1.summary[[2]][c('2.5%', '50%', '97.5%')]
post_means<-colMeans(posterior_samples.all)

ve.simple<- round(c(post_means, ci[,1]), 2) 
return(ve.simple)
}