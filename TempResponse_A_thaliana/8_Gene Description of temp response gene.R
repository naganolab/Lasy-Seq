##############
### script ###
##############


# Addition of Gene Description information  (For Table S2) ---------

load("R/Araport11_genes201606transcript_rep_ERCC_Virus7457_desc")

#Genes correlated to the temperature on sampling day
d<-read.csv("R/GO/corcoef005genes_d0_2921.csv")
gene<-as.character(d[,1])
des2<-des[gene,]
res<-cbind(d, des2)
write.csv(res, "R/GO/des_q01c005_d0.csv")


#Genes correlated to the temperature on 1day prior to sampling
d<-read.csv("R/GO/corcoef005genes_d1_435.csv")
gene<-as.character(d[,1])
des2<-des[gene,]
res<-cbind(d, des2)
write.csv(res, "R/GO/des_q01c005_d1.csv")


#Genes correlated to the temperature on 2day prior to sampling
d<-read.csv("R/GO/corcoef005genes_d2_371.csv")
gene<-as.character(d[,1])
des2<-des[gene,]
res<-cbind(d, des2)
write.csv(res, "R/GO/des_q01c005_d2.csv")


#Genes correlated to the temperature on 3day prior to sampling
d<-read.csv("R/GO/corcoef005genes_d3_8.csv")
gene<-as.character(d[,1])
des2<-des[gene,]
res<-cbind(d, des2)
write.csv(res, "des_q01c005_d3.csv")

