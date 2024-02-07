#Script used to generate the file time_comp.Rdata relative to Figure 2 of the paper 
#Comparing and Updating R packages using MCMC Algorithms for Linear Inverse Modeling of Metabolic Networks 

library(limSolve)
library(samplelim)
library(tidyverse)

timecomp<-data.frame(modele=character(),nUnknowns=integer(),sample_size=integer(),func=character(),sample_index=integer(),time=double())

lim<-df2lim("DeclarationFileBOWF-short.txt")
modele="BOWF-short"
nUnknowns=lim$NUnknowns

sample_sizes<-c(50,100,500,1000,5000,10000,50000)

row_number=1

for (n in ncomp){
  for (func in c("xsample()","rlim()")){
    for (i in 1:10){
      print(paste(n,func,i))
      timecomp[row_number,1]=modele
      timecomp[row_number,2]=nUnknowns
      timecomp[row_number,3]=n
      timecomp[row_number,4]=func
      timecomp[row_number,5]=i
      
      if (func=="xsample"){
        timecomp[row_number,6]=system.time(xsample(E=lim$A,F=lim$B,G=lim$G,H=lim$H, iter=n))[[3]]
      }else{
        timecomp[row_number,6]=system.time(rlim(lim,nsamp=n))[[3]]
      }
      row_number<-row_number+1
      
    } 
    
    
    
  }
    
    
}

# save(timecomp,file="time_comp.Rdata")

summarised_timecomp<-timecomp %>% group_by(sample_size,func) %>% mutate(time=mean(time))%>%subset(select=c("sample_size","func","time"))%>%distinct()

#Plot the corresponding graph
p<-ggplot(data=summarised_timecomp,aes(x=sample_size,y=time,group=func))+geom_line(aes(linetype=func,color=func)) +xlab("n")+ylab("Computation time (s)")
p+scale_color_brewer(palette="Dark2")+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),legend.text = element_text(size=14))

