##############
### script ###
##############

# temperature and expression -------------------------------------------------

# read rpm and description files
load("R/180509-3_rpm_Ath")   #rpm
load("R/Araport11_genes201606transcript_rep_ERCC_Virus7457_desc")   #gene description

#data read
rpmori<-rpm
rpm2<-rpmori[des$NormalizationGroup=="data",]
rpm<-log10(rpm2+1)

#temperature setting file
d<-read.csv("R/Temperature/Sample_temp.csv")

#temperature for each day
t0<-d$Temp_sample      #temperature on the sampling day
t1<-d$Temp_pre1day     #temperature on the 1 day prior to sampling
t2<-d$Temp_pre2day     #temperature on the 2 days prior to sampling
t3<-d$Temp_pre3day     #temperature on the 3 days prior to sampling


#plot ----------------------------------------------------------
pdf(file="R/TempSetting/tempcorrelation.pdf")
plot(t0,t1, main="current_pre1d temp q=0.049")
plot(t0,t2, main="current_pre2d temp q=0.00065")
plot(t0,t3, main="current_pre3d temp q=0.19")
dev.off()

#Cor.test----------------------------------------------------------------
cor.test(t0,t1) 
cor.out=unlist(cor.test(t0,t1))
as.matrix(cor.out)
out1<-cor.out

cor.test(t0,t2)
cor.out=unlist(cor.test(t0,t2))
as.matrix(cor.out)
out2<-cor.out

cor.test(t0,t3) 
cor.out=unlist(cor.test(t0,t3))
as.matrix(cor.out)
out3<-cor.out

out<-rbind(out1, out2, out3)
out<-as.matrix(out)
p.value<-out[,3] #p.value
q.value<-p.adjust(p.value, method="BH")
q.value<-as.matrix(q.value)

out<-cbind(out, q.value)
  
write.csv(out,"R/TempSetting/tempcorrelation.csv")

#Plot of temperature settings --------------------------------------------
d<-read.csv("R/TempSetting/TemperatureSettings.csv")

#Temperature on day8-21 in 3 incubators
ExpA<-d[16:43,2]        #incubator 1 (temperature setting1)
ExpB<-d[16:43,3]        #incubator 2 (temperature setting2)
ExpC<-d[16:43,4]        #incubator 3 (temperature setting3)

pdf("R/TempSetting/temperature settings.pdf")
plot(ExpA, type="l", lwd=2, col="#87ceeb")
points(ExpB, type="l", lty=("dashed"), lwd=2, col="#ff69b4")
points(ExpC, type="l", lty=("dotted"), lwd=2, col="#808000")  

dev.off()





