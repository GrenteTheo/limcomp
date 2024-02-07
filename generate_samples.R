
samples<-data.frame(modele=character(),nUnknowns=integer(),sampleSize=integer(),algo=character(),sample_index=integer(),name=character())
sample_size=20000
rowNumber<-1

lim<-df2lim("DeclarationFileBOWF-short.txt")
modele="BOWF-short"
nUnknowns=lim$NUnknowns
folder="samples/"
for (i in 1:10){
  for (alg in c("rlim()","xsample()")){
    samples[rowNumber,1]=modele
    samples[rowNumber,2]=nUnknowns
    samples[rowNumber,3]=sample_size
    samples[rowNumber,4]=alg
    samples[rowNumber,5]=i
    
    
    print(paste(alg,i))
    name<-paste(folder,modele,"_",alg,"_",sample_size,"_",i,".csv",sep="")
    
    if (alg=="xsample()"){
      raw_data<-xsample(E=lim$A,F=lim$B,G=lim$G,H=lim$H, iter=sample_size)
    }else{
      raw_data<-rlim(lim,nsamp=sample_size)
    }
    write.csv(raw_data,name)
    samples[rowNumber,6]=name
    rowNumber<-rowNumber+1
 
  }
}

save(samples,file = "samples_index.Rdata")
