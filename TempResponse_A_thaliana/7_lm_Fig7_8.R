##############
### script ###
##############

# temperature and expression -------------------------------------------------

# rpm and gene description files
load("R/190322-3_rpm_Ath")   #rpm
load("R/Araport11_genes201606transcript_rep_ERCC_Virus7457_desc")   #des

#data read
rpmori<-rpm
rpm<-rpmori[des$NormalizationGroup=="data",]
rpm <- rpm[rowSums(rpm[,types==0]>0)>0, ]
log.rpm <- log10(rpm+1)
rcol<-log.rpm[,d$type=="Col0"]    #samples for Lasy-Seq

#Temprature settings
d<-read.csv("R/Temperature/Sample_temp.csv")
temp<-d$Temp_sample
types <- as.numeric(d$type) - 1
dcol<-subset(d, d$type=="Col0")  #samples for Lasy-Seq
d0<-dcol$Temp_sample        #temperature on the sampling day
d1<-dcol$Temp_pre1day       #temperature on the 1 day prior to sampling
d2<-dcol$Temp_pre2day       #temperature on the 2 days prior to sampling
d3<-dcol$Temp_pre3day       #temperature on the 3 days prior to sampling


#Genes correlated to the temperature on sampling day--------------------------

pdf(file="R/Temperature/day0_q01genes.pdf")
genes<-read.csv("R/Temperature/corcoef005genes_d02921.csv")
gene<-genes[,1]
g=as.character(gene)
par(oma = c(0, 0, 0, 0))  
par(mfrow=c(2,2)) 
i=1
for(i in 1:length(g)){
  j=g[i]
  j2=  subset(des[,2], rownames(des)==j)
  plot(d0, rcol[j,], main=paste(j2,"_",j),cex.main=0.8, xlab="sampling temperature") 
  lm.obj<-lm(rcol[j,]~d0)
  abline(lm.obj,col=2)
  mtext(substitute(paste(p-value,"=",x),list(x=summary(lm.obj)$coefficients[8],digits=5)), side = 3,col=2,adj=0)
  
  plot(d1, rcol[j,], main=paste(j2,"_",j),cex.main=0.8, xlab="pre 1 day temperature") 
  lm.obj<-lm(rcol[j,]~d1)
  abline(lm.obj,col=2)
  mtext(substitute(paste(p-value,"=",x),list(x=summary(lm.obj)$coefficients[8],digits=5)), side = 3,col=2,adj=0)
  
  plot(d2, rcol[j,], main=paste(j2,"_",j),cex.main=0.8, xlab="pre 2 day temperature") 
  lm.obj<-lm(rcol[j,]~d2)
  abline(lm.obj,col=2)
  mtext(substitute(paste(p-value,"=",x),list(x=summary(lm.obj)$coefficients[8],digits=5)), side = 3,col=2,adj=0)
  
  plot(d3, rcol[j,], main=paste(j2,"_",j),cex.main=0.8, xlab="pre 3 day temperature") 
  lm.obj<-lm(rcol[j,]~d3)
  abline(lm.obj,col=2)
  mtext(substitute(paste(p-value,"=",x),list(x=summary(lm.obj)$coefficients[8],digits=5)), side = 3,col=2,adj=0)
  }
dev.off()

#Genes correlated to the temperature on 1day prior to sampling--------------------------

pdf(file="R/Temperature/pre1day_q01genes.pdf")
genes<-read.csv("R/Temperature/corcoef005genes_d1_435.csv")
gene<-genes[,1]
g=as.character(gene)
par(oma = c(0, 0, 0, 0))   
par(mfrow=c(2,2)) 
i=1
for(i in 1:length(g)){
  j=g[i]
  j2=  subset(des[,2], rownames(des)==j)
  plot(d0, rcol[j,], main=paste(j2,"_",j),cex.main=0.8, xlab="sampling temperature") 
  lm.obj<-lm(rcol[j,]~d0)
  abline(lm.obj,col=2)
  mtext(substitute(paste(p-value,"=",x),list(x=summary(lm.obj)$coefficients[8],digits=5)), side = 3,col=2,adj=0)
  
  plot(d1, rcol[j,], main=paste(j2,"_",j),cex.main=0.8, xlab="pre 1 day temperature") 
  lm.obj<-lm(rcol[j,]~d1)
  abline(lm.obj,col=2)
  mtext(substitute(paste(p-value,"=",x),list(x=summary(lm.obj)$coefficients[8],digits=5)), side = 3,col=2,adj=0)
  
  plot(d2, rcol[j,], main=paste(j2,"_",j),cex.main=0.8, xlab="pre 2 day temperature") 
  lm.obj<-lm(rcol[j,]~d2)
  abline(lm.obj,col=2)
  mtext(substitute(paste(p-value,"=",x),list(x=summary(lm.obj)$coefficients[8],digits=5)), side = 3,col=2,adj=0)
  
  plot(d3, rcol[j,], main=paste(j2,"_",j),cex.main=0.8, xlab="pre 3 day temperature") 
  lm.obj<-lm(rcol[j,]~d3)
  abline(lm.obj,col=2)
  mtext(substitute(paste(p-value,"=",x),list(x=summary(lm.obj)$coefficients[8],digits=5)), side = 3,col=2,adj=0)
}
dev.off()



#Genes correlated to the temperature on 2day prior to sampling--------------------------

pdf(file="R/Temperature/pre2day_q01genes.pdf")
genes<-read.csv("R/Temperature/corcoef005genes_d2_371.csv")
gene<-genes[,1]
g=as.character(gene)
par(oma = c(0, 0, 0, 0))   
par(mfrow=c(2,2)) 
i=1
for(i in 1:length(g)){
  j=g[i]
  j2=  subset(des[,2], rownames(des)==j)
  plot(d0, rcol[j,], main=paste(j2,"_",j),cex.main=0.8, xlab="sampling temperature") 
  lm.obj<-lm(rcol[j,]~d0)
  abline(lm.obj,col=2)
  mtext(substitute(paste(p-value,"=",x),list(x=summary(lm.obj)$coefficients[8],digits=5)), side = 3,col=2,adj=0)
  
  plot(d1, rcol[j,], main=paste(j2,"_",j),cex.main=0.8, xlab="pre 1 day temperature") 
  lm.obj<-lm(rcol[j,]~d1)
  abline(lm.obj,col=2)
  mtext(substitute(paste(p-value,"=",x),list(x=summary(lm.obj)$coefficients[8],digits=5)), side = 3,col=2,adj=0)
  
  plot(d2, rcol[j,], main=paste(j2,"_",j),cex.main=0.8, xlab="pre 2 day temperature") 
  lm.obj<-lm(rcol[j,]~d2)
  abline(lm.obj,col=2)
  mtext(substitute(paste(p-value,"=",x),list(x=summary(lm.obj)$coefficients[8],digits=5)), side = 3,col=2,adj=0)
  
  plot(d3, rcol[j,], main=paste(j2,"_",j),cex.main=0.8, xlab="pre 3 day temperature") 
  lm.obj<-lm(rcol[j,]~d3)
  abline(lm.obj,col=2)
  mtext(substitute(paste(p-value,"=",x),list(x=summary(lm.obj)$coefficients[8],digits=5)), side = 3,col=2,adj=0)
}
dev.off()


#Genes correlated to the temperature on 3day prior to sampling--------------------------

pdf(file="R/Temperature/pre3day_q01genes.pdf")
genes<-read.csv("R/Temperature/corcoef005genes_d3_8.csv") #all genes
gene<-genes[,1]
g=as.character(gene)
par(oma = c(0, 0, 0, 0))  
par(mfrow=c(2,2)) 
i=1
for(i in 1:length(g)){
  j=g[i]
  j2=  subset(des[,2], rownames(des)==j)
  plot(d0, rcol[j,], main=paste(j2,"_",j),cex.main=0.8, xlab="sampling temperature") 
  lm.obj<-lm(rcol[j,]~d0)
  abline(lm.obj,col=2)
  mtext(substitute(paste(p-value,"=",x),list(x=summary(lm.obj)$coefficients[8],digits=5)), side = 3,col=2,adj=0)
  
  plot(d1, rcol[j,], main=paste(j2,"_",j),cex.main=0.8, xlab="pre 1 day temperature") 
  lm.obj<-lm(rcol[j,]~d1)
  abline(lm.obj,col=2)
  mtext(substitute(paste(p-value,"=",x),list(x=summary(lm.obj)$coefficients[8],digits=5)), side = 3,col=2,adj=0)
  
  plot(d2, rcol[j,], main=paste(j2,"_",j),cex.main=0.8, xlab="pre 2 day temperature") 
  lm.obj<-lm(rcol[j,]~d2)
  abline(lm.obj,col=2)
  mtext(substitute(paste(p-value,"=",x),list(x=summary(lm.obj)$coefficients[8],digits=5)), side = 3,col=2,adj=0)
  
  plot(d3, rcol[j,], main=paste(j2,"_",j),cex.main=0.8, xlab="pre 3 day temperature") 
  lm.obj<-lm(rcol[j,]~d3)
  abline(lm.obj,col=2)
  mtext(substitute(paste(p-value,"=",x),list(x=summary(lm.obj)$coefficients[8],digits=5)), side = 3,col=2,adj=0)
}
dev.off()

