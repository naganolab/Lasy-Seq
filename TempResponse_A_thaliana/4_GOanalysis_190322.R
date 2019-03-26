#################
# preparation   #
#################

source("https://bioconductor.org/biocLite.R")
biocLite("GO.db")
library("GO.db")

# load the GO analysis functions
fn <- sprintf("R/GO/1_GOanalysis_functions_180613.R")
source(fn)

# load the table of Gene ID and GO (as ulg)
fn <- sprintf("R/GO/ulg.TAIR_170801") # for A. thaliana
load(fn) 


#################
#   Analysis    #
#################

# load an gene list (as g0~g3)
g0<-read.csv("R/GO/corcoef005genes_d0_2921.csv")
g1<-read.csv("R/GO/corcoef005genes_d1_435.csv")
g2<-read.csv("R/GO/corcoef005genes_d2_371.csv")
g3<-read.csv("R/GO/corcoef005genes_d3_8.csv")


#For genes which had significant correlation to the sampling day
gl<-g0[,1]
#write.csv(gl, "genelist.csv")
# fisher's exact test for all GO
cgt=ulg
gn.t=gl
gn.test=substr(gn.t, 1,9) 
result <- ng.mft(cgt, gn.test)   
ng.prepGOtestOutTable(result)
# output results as a csv file
write.csv(ng.prepGOtestOutTable(result), file="GO_q01cor005_d3.csv")
write.csv(result, file="Matrix_q01cor005_d3.csv")


#For genes which had significant correlation to 1-day prior to sampling
gl<-g1[,1]  #change for each date
cgt=ulg
gn.t=gl
gn.test=substr(gn.t, 1,9) 
result <- ng.mft(cgt, gn.test)   
ng.prepGOtestOutTable(result)
write.csv(ng.prepGOtestOutTable(result), file="GO_q01cor005_d1.csv")
write.csv(result, file="Matrix_q01cor005_d1.csv")


#genes which had significant correlation to 2-days prior to sampling
gl<-g2[,1]
cgt=ulg
gn.t=gl
gn.test=substr(gn.t, 1,9) 
result <- ng.mft(cgt, gn.test)   
ng.prepGOtestOutTable(result)
write.csv(ng.prepGOtestOutTable(result), file="GO_q01cor005_d2.csv")
write.csv(result, file="Matrix_q01cor005_d2.csv")


#genes which had significant correlation to 3-days prior to sampling
gl<-g3[,1]
cgt=ulg
gn.t=gl
gn.test=substr(gn.t, 1,9) 
result <- ng.mft(cgt, gn.test)   
ng.prepGOtestOutTable(result)
write.csv(ng.prepGOtestOutTable(result), file="GO_q01cor005_d3.csv")
write.csv(result, file="Matrix_q01cor005_d3.csv")





##############
#   Result   #
##############

#d0
Adjusted P value       ID           Description                                A & B  A      B       U      
GO:0005622 "1.39035316855299e-08" "GO:0005622" "intracellular"                            "1594" "2145" "20460" "30260"
GO:0005623 "5.24383532764452e-10" "GO:0005623" "cell"                                     "1730" "2145" "22405" "30260"
GO:0043231 "0.00690658748422726"  "GO:0043231" "intracellular membrane-bounded organelle" "1377" "2145" "17937" "30260"
GO:0043226 "0.0198876130497534"   "GO:0043226" "organelle"                                "1393" "2145" "18235" "30260"
GO:0043227 "0.00760345156146033"  "GO:0043227" "membrane-bounded organelle"               "1377" "2145" "17943" "30260"
GO:0043229 "0.0184081715890685"   "GO:0043229" "intracellular organelle"                  "1393" "2145" "18230" "30260"
GO:0044424 "1.89055942592773e-08" "GO:0044424" "intracellular part"                       "1591" "2145" "20428" "30260"
GO:0044464 "4.4743826235738e-10"  "GO:0044464" "cell part"                                "1730" "2145" "22399" "30260"

#d1
Adjusted P value      ID           Description          A & B A     B       U      
GO:0005622 "0.00422223076084171" "GO:0005622" "intracellular"      "251" "313" "20460" "30260"
GO:0005623 "0.0189726444990469"  "GO:0005623" "cell"               "266" "313" "22405" "30260"
GO:0044424 "0.00682058476216159" "GO:0044424" "intracellular part" "250" "313" "20428" "30260"
GO:0044464 "0.0182922776139558"  "GO:0044464" "cell part"          "266" "313" "22399" "30260"

#d2
Adjusted P value ID Description A & B A B U

#d3
Adjusted P value ID Description A & B A B U




