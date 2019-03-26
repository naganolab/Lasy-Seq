
##############
### script ###
##############

# set parameters -----------------------------------------------------
dir.input <- "R"
dir.output <- "R"
dir.rawdata <- "result"

fn.description <- "141103-1_GeneDescription_Osa"
fn.attribute <- "SampleAttribute_Osa.txt"

exec.date <- "190323"


fn.rawcnt <- sprintf("%s-2_rawcnt_Osa", exec.date)
fn.rpm <- sprintf("%s-3_rpm_Osa", exec.date)
fn.cvrd1 <- sprintf("%s-4_cvrd1_Osa", exec.date)
fn.cvrd3 <- sprintf("%s-5_cvrd3_Osa", exec.date)


#load sample attribute --------------
fn <- sprintf("%s/%s", dir.input, fn.attribute)
at <- read.delim(fn, header=T, as.is=T)
#rownames(at) <- sprintf(at[,"name"], at[,"library"])
rownames(at) <- sprintf("%s_%03d", at[,"library"], at[,"index"])



#load Ahg transcript description --------------------
fn <- sprintf("%s/%s", dir.input, fn.description)
load(fn)

#set data read row --------------
drr <- des[,"NormalizationGroup"]=="data"
names(drr) <- rownames(des)

# save/load expression data table-------------------------------------
#rpkm <- NULL
rawcnt <- matrix(0, nrow=nrow(des), ncol=nrow(at))
cvrd1 <- rawcnt
cvrd3 <- rawcnt

for(i in 1:nrow(at)){
  #for(i in 580:590){
  fn <- sprintf("%s/%s/%s/Index%03d.genes.results",
                dir.output, dir.rawdata, at[i,"library"], at[i,"index"])
  tmp <- read.delim(fn, as.is=T, comment.char="#")
  rawcnt[,i] <- tmp[,"expected_count"]
  if(is.element("coverage.depth..1.",colnames(tmp))){
    cvrd1[,i] <- tmp[,"coverage.depth..1."]
    cvrd3[,i] <- tmp[,"coverage.depth..3."]
  } 
  
  cat(i)
  cat("\n")
}

colnames(rawcnt) <- rownames(at)
colnames(cvrd1) <- rownames(at)
colnames(cvrd3) <- rownames(at)
rownames(rawcnt) <- rownames(des)
rownames(cvrd1) <- rownames(des)
rownames(cvrd3) <- rownames(des)

#calc. rpm
rce <- rawcnt[des[,"NormalizationGroup"]=="data", ]
tmp <- colSums(rce)/10^6
rce2 <- t(t(rce)/tmp)

rcv <- rawcnt[des[,"NormalizationGroup"]=="virus", ]
tmp <- colSums(rce)/10^6
rcv2 <- t(t(rcv)/tmp)

rcb <- rawcnt[des[,"NormalizationGroup"]=="rRNA", ]
tmp <- colSums(rcb)/10^6
rcb2 <- t(t(rcb)/tmp)

rcc <- rawcnt[des[,"NormalizationGroup"]=="ercc", ]
tmp <- colSums(rcc)/10^6
rcc2 <- t(t(rcc)/tmp)

rpm <- rbind(rce2, rcv2, rcb2, rcc2)
rpm[is.nan(rpm)] <- 0
rpm <- rpm[rownames(des),]

# subtruction of indexing contamination
ulib <- unique(at[,"library"])

rpmtmp <- NULL
for(i in ulib){
  extmp <- rpm[, at[,"library"]==i]
  tmp <- rowSums(extmp)*0.0005
  extmp <- extmp-tmp
  rpmtmp <- cbind(rpmtmp, extmp)
}
rpmtmp[rpmtmp<0] <- 0

rpm <- rpmtmp[,rownames(at)]
cvrd1[rpm==0] <- 0
cvrd3[rpm==0] <- 0

# save objects
save(rawcnt, file=sprintf("%s/%s", dir.output, fn.rawcnt))
save(rpm, file=sprintf("%s/%s", dir.output, fn.rpm))
save(cvrd1, file=sprintf("%s/%s", dir.output, fn.cvrd1))
save(cvrd3, file=sprintf("%s/%s", dir.output, fn.cvrd3))







