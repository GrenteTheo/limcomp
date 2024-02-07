#Script used to generate the file Example2_diag.Rdata relative to Table 2 of the paper 
#Comparing and Updating R packages using MCMC Algorithms for Linear Inverse Modeling of Metabolic Networks 

library(coda)
library(statip)
library(pracma)
library(stats)
library(tidyverse)

load(file="samples_index.Rdata")

lim<-df2lim("DeclarationFileBOWF-short.txt")

theorical_ranges<-lim.ranges(lim)


current_row<-list(modele=character(),nUnknowns=integer(),sampleSize=integer(),algo=character(),jmp_type=character(),jmp=double(),sampleIndex=integer(),fluxIndex=double(),flux_min=double(),flux_max=double(),flux_mean=double(),flux_median=double(),flux_Q1=double(),flux_Q3=double(),RL_M=double(),RL_N=double(),RL_Nmin=double(),RL_I=double(),RL_k=double(),Geweke=double(),HD=double(),ESS=double(),sample=list(),theorical_min=double(),theorical_max=double(),rangeCoverage=double())

diag_dataset=list()

for (n in 1:nrow(samples)){
  print(n)
  row_smp<-samples[n,]
  current_row["modele"]=row_smp[1]
  nUnknowns=row_smp[2]
  current_row["nUnknowns"]=nUnknowns
  current_row["sampleSize"]=row_smp[3]
  alg=row_smp[4]
  current_row["algo"]=alg
  i=row_bench[5]
  current_row["sampleIndex"]=i
  sample<-read.csv(row_bench[[6]],header = TRUE)[,2:(nUnknowns+1)]
  for (j in 1:nUnknowns){
    current_row["fluxIndex"]=j
    flux_sample<-c(sample[,j])
    current_row["sample"]<-list(flux_sample)
    sample_quantiles<-quantile(flux_sample,names=FALSE)
    current_row["flux_min"]<-sample_quantiles[1]
    current_row["flux_Q1"]<-sample_quantiles[2]
    current_row["flux_median"]<-sample_quantiles[3]
    current_row["flux_Q3"]<-sample_quantiles[4]
    current_row["flux_max"]<-sample_quantiles[5]
    current_row["flux_mean"]<-mean(flux_sample)
    
    #Raftery Lewis
      
    data<-mcmc(data= flux_sample, start = 1, end = niter, thin = 1)
    raftery_lewis<-raftery.diag(data, q=0.025, r=0.005, s=0.95, converge.eps=0.001)
    res_rl<-raftery_lewis[["resmatrix"]]
    current_row["RL_M"]<-res_rl[[1,1]]
    current_row["RL_N"]<-res_rl[[1,2]]
    current_row["RL_Nmin"]<-res_rl[[1,3]]
    current_row["RL_I"]<-res_rl[[1,4]]
    current_row["RL_k"]<-res_rl[[1,5]]
    
    #Geweke
    
    current_row["Geweke"]<-geweke.diag(data, frac1=0.1, frac2=0.5)$z
    
    #Hellinger Distance
    r1<-as.integer(niter/3)
    r2<-as.integer(2*niter/3)
    p1<-as.vector(flux_sample[1:r1])
    p1_dens<-densityfun(p1)
    p3<-as.vector(flux_sample[r2:niter])
    p3_dens<-densityfun(p3)
    fun <-function(x) (sqrt(p1_dens(x))-sqrt(p3_dens(x)))**2
    current_row["HD"]<-sqrt(0.5*integral(fun=fun,-Inf,Inf) )
    
    #ESS
    current_row["ESS"]<-effectiveSize(data)
    
    #Range Coverage
    minimum=theorical_ranges[[j,1]]
    current_row["theorical_min"]=minimum
    maximum=theorical_ranges[[j,2]]
    current_row["theorical_max"]=maximum
    if (maximum==minimum){
      current_row["rangeCoverage"]<-NA
    }else{
      current_row["rangeCoverage"]<-(maximum-minimum)/(current_row["flux_max"]-current_row["flux_min"])
      }
    
    # Add to big dataset
    diag_dataset<-rbind(diag_dataset,current_row)
  }
  
  }


