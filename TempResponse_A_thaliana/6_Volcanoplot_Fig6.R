##############
### script ###
##############

# temperature and expression -------------------------------------------------

# read rpm and description files
load("R/190322-3_rpm_Ath")   #rpm
load("R/Araport11_genes201606transcript_rep_ERCC_Virus7457_desc")   #des

# data read
rpmori<-rpm
rpm<-rpmori[des$NormalizationGroup=="data",]
rpm <- rpm[rowSums(rpm[,types==0]>0)>0, ]
log.rpm <- log10(rpm+1)
#samples for Lasy-Seq
rcol<-log.rpm[,d$type=="Col0"]

# tempearature setting of the sample
d<-read.csv("R/Temperature/Sample_temp.csv")
temp<-d$Temp_sample
types <- as.numeric(d$type) - 1
#samples for Lasy-Seq
dcol<-subset(d, d$type=="Col0")
d0<-dcol$Temp_sample        #temperature on the sampling day
d1<-dcol$Temp_pre1day       #temperature on the 1 day prior to sampling
d2<-dcol$Temp_pre2day       #temperature on the 2 days prior to sampling
d3<-dcol$Temp_pre3day       #temperature on the 3 days prior to sampling


# lm ----------------------------------------------------------------------------------------

e.out<-matrix(0, nrow=5, ncol=nrow(rcol))
colnames(e.out)<-rownames(rcol)
rownames(e.out)<-c("Intercept_Estimate","d0_Estimate","d1_Estimate","d2_Estimate","d3_Estimate")
p.out<-matrix(0, nrow=5, ncol=nrow(rcol))
colnames(p.out)<-rownames(rcol)
rownames(p.out)<-c("Intercept_pvalue","d0_pvalue","d1_pvalue","d2_pvalue","d3_pvalue")

r.out<-matrix(0, nrow=1, ncol=nrow(rcol)) 
colnames(r.out)<-rownames(rcol)
rownames(r.out)<-c("r.squared of the full model")

for(i in 1:nrow(rcol)){
l<-lm(rcol[i,] ~ scale(d0)+scale(d1)+scale(d2)+scale(d3))
s<-summary(l)
Estimate<-s$coefficients[,1]
p.value<-s$coefficients[,4]
R.squared<-s$r.squared
e.out[,i]<-Estimate
p.out[,i]<-p.value
r.out[,i]<-R.squared
}

# FDR control 
q_Intercept<-p.adjust(p.out[1,], "fdr")
q_d0<-p.adjust(p.out[2,], "fdr")
q_d1<-p.adjust(p.out[3,], "fdr")
q_d2<-p.adjust(p.out[4,], "fdr")
q_d3<-p.adjust(p.out[5,], "fdr")

# save
write.csv(q_d0, "R/Temperature/q_d0.csv")
write.csv(q_d1, "R/Temperature/q_d1.csv")
write.csv(q_d2, "R/Temperature/q_d2.csv")
write.csv(q_d3, "R/Temperature/q_d3.csv")

# number of genes significuntly correlated to the temperature  (q<0.1, FDR)
q0s<-subset(q_d0, q_d0<0.1)  #temperature on the 1 day prior to sampling 
length(q0s)                  #3007
q1s<-subset(q_d1, q_d1<0.1)  #temperature on the 1 day prior to sampling
length(q1s)                  #456
q2s<-subset(q_d2, q_d2<0.1)  #temperature on the 2 days prior to sampling
length(q2s)                  #376
q3s<-subset(q_d3, q_d3<0.1)  #temperature on the 3 days prior to sampling
length(q3s)                  #8


#histogram for Figure 6-----------------------------------------------------------------------

pdf(file="R/Temperature/e.out_hist.pdf")
par(oma = c(0, 0, 0, 0))   # 下・左・上・右の順で余白を設定
par(mfrow=c(2,2)) #画面を 2 × 2 に分割
hist(e.out[2,], xlim=c(-1,1),breaks = 40, ylim=c(0,9000)) #cor coef for sampling day
hist(e.out[3,], xlim=c(-1,1), ylim=c(0,9000))  #pre-1 day 
hist(e.out[4,], xlim=c(-1,1), ylim=c(0,9000))  #pre-2 day 
hist(e.out[5,], xlim=c(-1,1), ylim=c(0,9000))  #pre-3 day 
dev.off()


#volcano plot with p-value for Figure 6-------------------------------------------------------

pdf(file="R/Temperature/Volcano_FDR005orange_lm.pdf")

par(oma = c(0, 0, 0, 0))  
par(mfrow=c(2,2)) 

# Volcano plot for temperature on the sampling day 
y<--log10(p.out[2,])
x<-e.out[2,]
plot(x,y, xlim=c(-1,1),col=rgb(0, 0, 0, alpha=0.4), ylim=c(0,18),xlab="correlation coefficients", ylab="-log10(q-value)", main="sampling day")
y2<-y[q_d0<0.1]
x2<-x[q_d0<0.1]
points(x2,y2, col = "red") #significant genes
points(x2[x2>-0.05 & x2< 0.05],y2[x2>-0.05 & x2< 0.05], col = "Gold") #FDR(q)<0.1 and correlation coefficients<0.3
length(x2[x2>-0.05 & x2< 0.05]) 
a<-x2[x2> 0.05]
b<-x2[x2< -0.05] 
a2<-as.matrix(a)
b2<-as.matrix(b)
dim(rbind(a2,b2)) #2921, 1   
write.csv(rbind(a2,b2), "R/Temperature/corcoef005genes_d02921.csv")

# Volcano plot for temperature on the 1 day prior to sampling 
y<--log10(p.out[3,])
x<-e.out[3,]
plot(x,y, xlim=c(-1,1),col=rgb(0, 0, 0, alpha=0.4), ylim=c(0,18),xlab="correlation coefficients", ylab="-log10(q-value)", main="pre-1 day")
y2<-y[q_d1<0.1]
x2<-x[q_d1<0.1]
points(x2,y2, col = "red")
points(x2[x2>-0.05 & x2< 0.05],y2[x2>-0.05 & x2< 0.05], col = "Gold")
length(x2[x2>-0.05 & x2< 0.05])
a<-x2[x2> 0.05] 
b<-x2[x2< -0.05] 
a2<-as.matrix(a)
b2<-as.matrix(b)
dim(rbind(a2,b2)) #435, 1   
write.csv(rbind(a2,b2), "R/Temperature/corcoef005genes_d1_435.csv")

# Volcano plot for temperature on the 2 days prior to sampling 
y<--log10(p.out[4,])
x<-e.out[4,]
plot(x,y, xlim=c(-1,1),col=rgb(0, 0, 0, alpha=0.4), ylim=c(0,18),xlab="correlation coefficients", ylab="-log10(q-value)", main="pre-2 day")
y2<-y[q_d2<0.1]
x2<-x[q_d2<0.1]
points(x2,y2, col = "red")
points(x2[x2>-0.05 & x2< 0.05],y2[x2>-0.05 & x2< 0.05], col = "Gold")
length(x2[x2>-0.05 & x2< 0.05]) 
a<-x2[x2> 0.05] 
b<-x2[x2< -0.05] 
a2<-as.matrix(a)
b2<-as.matrix(b)
dim(rbind(a2,b2)) #371, 1  
write.csv(rbind(a2,b2), "R/Temperature/corcoef005genes_d2_371.csv")

# Volcano plot for temperature on the 3 days prior to sampling 
y<--log10(p.out[5,])
x<-e.out[5,]
plot(x,y, xlim=c(-1,1),col=rgb(0, 0, 0, alpha=0.4), ylim=c(0,18),xlab="correlation coefficients", ylab="-log10(q-value)", main="pre-3 day")
y2<-y[q_d3<0.1]
x2<-x[q_d3<0.1]
points(x2,y2, col = "red")
points(x2[x2>-0.05 & x2< 0.05],y2[x2>-0.05 & x2< 0.05], col = "Gold")
length(x2[x2>-0.05 & x2< 0.05]) 
a<-x2[x2> 0.05]
b<-x2[x2< -0.05] 
a2<-as.matrix(a)
b2<-as.matrix(b)
dim(rbind(a2,b2)) #8, 1   
write.csv(rbind(a2,b2), "R/Temperature/corcoef005genes_d3_8.csv")

dev.off()


