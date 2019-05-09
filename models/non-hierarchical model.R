non_hierarchical_func<-function(){
  model_string<-"
model{
  
  for(i in 1:2){
  
  #Likelihood
  y[i, 1:14] ~ dmulti(p[i, 1:14], m[i])
  
  #Multionmial Logistic Regression    
  for(j in 1:13){
  p_temp[i,j] <- exp(beta0[j] + beta1[j]*vax[i])
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
  for(j in 1:13){
  
  beta0[j] ~ dnorm(0, 0.0001)
  beta1[j] ~ dnorm(0, 0.0001)
  
  }
  
  #Overall Vaccine Effect
  overall_vax_effect <- 100*(1 - (1 - p[2,14])/(1 - p[1,14]))
  
}
"


##############################################################
#Model Fitting
##############################################################
jags.inits1 <- function(){
  list(".RNG.seed"=c(1), ".RNG.name"='base::Wichmann-Hill')
}

model_spec<-textConnection(model_string)
model_jags<-jags.model(model_spec, 
                       inits = jags.inits1,
                       data=list('y' = y,
                                 'm' = m,
                                 'vax' = vax),
                       n.adapt=10000, 
                       n.chains=3)

params<-c('overall_vax_effect',
          'sero_vax_effect')


##############################################
#Posterior Sampling
##############################################
posterior_samples<-coda.samples(model_jags, 
                                params, 
                                n.iter=100000)
posterior_samples.all<-do.call(rbind,posterior_samples)

###################################################
#POSTERIOR INFERENCE
###################################################
post_means<-colMeans(posterior_samples.all)
ci<-t(hdi(posterior_samples.all, credMass = 0.95))
ci<-round(ci,1)
post_means<-round(post_means,1)

st.labs<-as.character(unique(d1$st))

yrange<-range(ci)

overall.VE<-c(post_means[1], ci[1,])
st.VE<- cbind(post_means[-1], ci[-1,])


#install.packages('rmeta')
library(rmeta)
summary_data <- 
  structure(list(
    mean  = c(NA, NA,post_means[2:14],NA,post_means[1]), 
    lower = c(NA, NA,ci[2:14,1],NA,ci[1,1]),
    upper = c(NA, NA,ci[2:14,2],NA,ci[1,2]),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -11L), 
    class = "data.frame"))

if (outcome.var=='n_s11_pp'){
  ci['sero_vax_effect[8]',]<-NA
  post_means['sero_vax_effect[8]']<-NA
}
tabletext<-cbind(
  c("", "Serotype", st.labs, NA, "All"),
  c("Vaccine", "(N=42240)", d1[,outcome.var][d1$vax==1], NA, sum(d1[,outcome.var][d1$vax==1])),
  c("Control", "(N=42256)",  d1[,outcome.var][d1$vax==0], NA, sum(d1[,outcome.var][d1$vax==0])),
  c("", "VE", paste0(post_means[2:14],'% (', ci[2:14,1],'%, ',ci[2:14,2] ,'%' ,')'), 
    NA, paste0(post_means[1],'% (', ci[1,1],'%, ', ci[1,2],'%',')' ) )
)
    
  
res.list<-list('tabletext'=tabletext, 'summary_data'=summary_data,'st.VE'=st.VE)
return(res.list)

}