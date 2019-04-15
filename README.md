# capita-hierarchical
Re-analysis of CAPITA trial using hierarchical Bayesian model
This code re-analyzes the data in Bonten et al NEJM; table S11 (non-invasive, non-bacteremic pneumonia).

The capita1.R file is used to call 3 models:
A model with no serotype-specific effect (overall efficacy); a model with a serotype-specific effect but no pooling; 
and a hierarchical model where there is a hyperprior on the serotype-specific effect
