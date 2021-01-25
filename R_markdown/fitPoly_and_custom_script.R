# Script for going from AxiomGT1.summary.txt to compared dosage calls
# Author: Jeekin Lau
# Custom script at end used for comparing probes of Axiom WagRHSNP 68k array

library(SNPolisher)

fitTetra_Input(summaryFile="AxiomGT1.summary.txt",
               output.file="AxiomGT1.summary.fitTetra.txt")
			   
library(fitPoly)
library(doParallel)
data1<-read.delim("AxiomGT1.summary.fitTetra.txt",
                 sep="\t",stringsAsFactors=F)
				 
saveMarkerModels(4, markers=NA, data=data1, diplo=NULL, select=TRUE,
                 diploselect=TRUE, pop.parents=NULL, population=NULL, parentalPriors=NULL,
                 samplePriors=NULL, startmeans=NULL, maxiter=40, maxn.bin=200, nbin=200,
                 sd.threshold=0.1, p.threshold=0.99, call.threshold=0.6, peak.threshold=0.85,
                 try.HW=TRUE, dip.filter=1, sd.target=NA,
                 filePrefix=paste0(getwd(),"/plate13"), rdaFiles=F, allModelsFile=T,
                 plot ="none", plot.type="png", ncores=6)
				 
				 
				 
library(data.table)
calls<- as.matrix(fread("plate13_scores.dat", select = c(1:3,12)))
header_calls<-calls[1:1000,1:ncol(calls)]

num_ind<-96
ind<-calls[1:num_ind,3]
num_markers<-nrow(calls)/num_ind

markers<-matrix(,num_markers,1)
colnames(markers)<-"Probes_ID"

for (i in 1:nrow(markers)){
  markers[i,1]<-calls[i*num_ind,2]
  print(i)
}


genocalls<-matrix(, num_markers, num_ind+1)
colnames(genocalls)<-c("Probes_ID",ind)
genocalls[,1]<-markers

stringcall<-calls[,4]

genocalls<-t(genocalls)
genocalls[2:nrow(genocalls),1:ncol(genocalls)]<-stringcall
genocalls<-t(genocalls)


### order to right order
genocall_order<-as.matrix(read.csv("array_snps_flanking_order.csv"))

marker_col<-matrix(,nrow(genocalls),1)
genocalls2<-cbind(marker_col,genocalls)

for (a in 1:nrow(genocalls2)){
  probe<-genocalls2[a,2]
  genocalls2[a,1]<-genocall_order[which(genocall_order[,1]==probe),2]
  print(a)
}

genocalls3<-genocalls2[order(genocalls2[,1]),]

genocalls<-genocalls3


#compares the probes

compare_probes<-array(,dim=c(nrow(genocalls)/2,ncol(genocalls),2))

for (j in 1:nrow(compare_probes)){
  for (k in 1:ncol(compare_probes)){
    compare_probes[j,k,1]<-genocalls[j*2-1,k]
    compare_probes[j,k,2]<-genocalls[j*2,k]
    } 
  print(j)
}

write.csv(compare_probes[,,1],"probe1.csv", col.names=T, row.names = F)
write.csv(compare_probes[,,2],"probe2.csv", col.names=T, row.names = F)


probe1<-as.matrix(read.csv("probe1.csv", header=T))
probe2<-as.matrix(read.csv("probe2.csv", header=T))

compare_probes<-array(,dim=c(nrow(probe1),ncol(probe1),2))
compare_probes[,,1]<-probe1
compare_probes[,,2]<-probe2

compared_calls<-matrix(,nrow(compare_probes),ncol(compare_probes))

compare_probes[is.na(compare_probes)]<-9

for (l in 1:nrow(compared_calls)){
  for(m in 1:ncol(compared_calls)){
    p1<-compare_probes[l,m,1]
    p2<-compare_probes[l,m,2]
    
    ifelse(p1!=9 & p2!=9 & p1==p2, compared_calls[l,m]<-p1,
           ifelse(p1==9 & p2==9, compared_calls[l,m]<-NA,       
                  ifelse(p1==9 & p2!=9,compared_calls[l,m]<-p2,       
                         ifelse(p1!=9 & p2==9,compared_calls[l,m]<-p1,              
                                ifelse(p1!=9 & p2!=9 & p1 != p2, compared_calls[l,m]<-NA,                       
                                       NA)))))   
  } 
  print(l)}


for (b in 1:nrow(compared_calls)){
    compared_calls[b,1]<-genocalls3[b*2-1,1]
  print(b)
}

colnames(compared_calls)<-colnames(genocalls3)

write.csv(compared_calls, "compared_calls.csv", col.names = T, row.names = F)



compared_calls<-matrix(,nrow(compare_probes),ncol(compare_probes))

#S if the two probes are the same
#D if the calls are different
#O if the calls have one that is the same




for (l in 1:nrow(compared_calls)){
  for(m in 1:ncol(compared_calls)){
    p1<-compare_probes[l,m,1]
    p2<-compare_probes[l,m,2]
    
    ifelse(p1!=9 & p2!=9 & p1==p2, compared_calls[l,m]<-"S",
           ifelse(p1==9 & p2==9, compared_calls[l,m]<-NA,       
                  ifelse(p1==9 & p2!=9,compared_calls[l,m]<-"O",       
                         ifelse(p1!=9 & p2==9,compared_calls[l,m]<-"O",              
                                ifelse(p1!=9 & p2!=9 & p1 != p2, compared_calls[l,m]<-"D",                       
                                       NA)))))            
    
  } 
  print(l)}



write.csv(compared_calls, "compared_calls_kind_counts.csv")