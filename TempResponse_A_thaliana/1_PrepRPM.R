##############
### script ###
##############

# set parameters -----------------------------------------------------
dir.input <- "R"
dir.output <- "R"
dir.rawdata <- "RSEMout"

fn.description <- "Araport11_genes201606transcript_rep_ERCC_Virus7457_desc"
fn.sample.index <- "190322_Sample-IndexPrimer.txt"
fn.sample.index.merged <- "190322_Sample-IndexPrimer_jmj_merged.txt"

merge.sample <- TRUE
exec.date <- "190322"
species <- "Ath"

fn.rawcnt <- sprintf("%s-2_rawcnt_%s", exec.date, species)
fn.rpm <- sprintf("%s-3_rpm_%s", exec.date, species)
fn.cvrd1 <- sprintf("%s-4_cvrd1_%s", exec.date, species)
fn.cvrd3 <- sprintf("%s-5_cvrd3_%s", exec.date, species)
fn.readall <- sprintf("%s-6_readall_%s", exec.date, species)
fn.readfw <- sprintf("%s-7_readfw_%s", exec.date, species)
fn.fwpct <- sprintf("%s-8_fwpct_%s", exec.date, species)
fn.rcsummary <- sprintf("%s-1_ReadCountSummary_%s.txt", exec.date, species)


#load sample attribute --------------
fn <- sprintf("%s/%s", dir.input, fn.sample.index)
at <- read.delim(fn, header=T, as.is=T)
rownames(at) <- sprintf("%s_%03d", at[,"library"], at[,"index"])

#load transcript description --------------------
fn <- sprintf("%s/%s", dir.input, fn.description)
load(fn)

#set data read row --------------
drr <- des[,"NormalizationGroup"]=="data"
names(drr) <- rownames(des)


# save/load expression data table-------------------------------------
rawcnt <- matrix(0, nrow=nrow(des), ncol=nrow(at))
colnames(rawcnt) <- rownames(at)
rownames(rawcnt) <- rownames(des)
cvrd1 <- rawcnt
cvrd3 <- rawcnt
readall <- rawcnt
readfw <- rawcnt
fwpct <- rawcnt
readcount <- matrix(0, nrow=nrow(at), ncol=4)
colnames(readcount) <- c("Total reads", "Unique sequences", "Mapped reads", "Mapped rate")

for(i in 1:nrow(at)){
  #for(i in 580:590){
  fn <- sprintf("%s/%s/%s/Index%03d.genes.results",
                dir.input, dir.rawdata, at[i,"library"], at[i,"index"])
  
  tmp <- read.delim(fn, nrows=5, fill=T, header=F, as.is=T)
  readcount[i,] <- c(tmp[2,2], tmp[3,2], tmp[4,2], tmp[5,2])
  #readcount[i,] <- c(tmp[2,2], tmp[2,2], tmp[3,2], tmp[4,2])
  
  tmp <- read.delim(fn, as.is=T, comment.char="#")
  rawcnt[,i] <- tmp[,"expected_count"]
  if(is.element("coverage.depth..1.",colnames(tmp))){
    cvrd1[,i] <- tmp[,"coverage.depth..1."]
    cvrd3[,i] <- tmp[,"coverage.depth..3."]
  } 
  if(is.element("mapped.reads.all.",colnames(tmp))){
    readall[,i] <- tmp[,"mapped.reads.all."]
    readfw[,i] <- tmp[,"mapped.reads.fw."]
    fwpct[,i] <- tmp[,"fw..."]
  } 
  
  #for paired-end data
  if(is.element("mapped.read.pairs.all.",colnames(tmp))){ 
    readall[,i] <- tmp[,"mapped.read.pairs.all."]
    readfw[,i] <- tmp[,"mapped.reads.fw."]
    fwpct[,i] <- tmp[,"fw..."]
  } 
  
  cat(i)
  cat("\n")
}



# subtruction of indexing contamination
ulib <- unique(at[,"library"])

rawcnttmp <- NULL
for(i in ulib){
  extmp <- rawcnt[, at[,"library"]==i]
  tmp <- rowSums(extmp)*0.0005
  extmp <- extmp-tmp
  rawcnttmp <- cbind(rawcnttmp, extmp)
}
rawcnttmp[rawcnttmp<0] <- 0

rawcnt <- rawcnttmp[,rownames(at)]
cvrd1[rawcnt==0] <- 0
cvrd3[rawcnt==0] <- 0
readall[rawcnt==0] <- 0
readfw[rawcnt==0] <- 0
fwpct[rawcnt==0] <- 0



# marge multi-fastq from 1 sample to 1 coloumn
if(merge.sample){

  uni.sampleID <- unique(at[,"sampleID"])
  
  # at
  at.tmp <- NULL
  for(i in 1:length(uni.sampleID)){
    tmp <- at[at[,"sampleID"] == uni.sampleID[i], ]
    
    if(is.vector(tmp)){
      at.tmp <- rbind(at.tmp, tmp)
    } else {
      at.tmp <- rbind(at.tmp, tmp[1,])
    }
  }
  at.ori <- at
  at <- at.tmp
  write.table(at, file=sprintf("%s/%s", dir.input, fn.sample.index.merged), sep="\t")
  at <- at.ori
  
  # rawcnt
  tmp <- rawcnt
  tmp.fun <- function(x){
    x2 <- aggregate(x, by=list(at[,"sampleID"]), FUN=sum)
    x3 <- x2[,2]
    names(x3) <- x2[,1]
    return(x3)
  }
  tmp2 <- apply(tmp, 1, FUN=tmp.fun)
  tmp3 <- t(tmp2)
  tmp4 <- tmp3[,uni.sampleID]
  rawcnt <- tmp4
}



#calc. rpm
rce <- rawcnt[des[,"NormalizationGroup"]=="data", ]
tmp <- colSums(rce)/10^6
rce2 <- t(t(rce)/tmp)
if(sum(des[,"NormalizationGroup"]=="data")==0){rce2<-NULL}

rcv <- rawcnt[des[,"NormalizationGroup"]=="virus", ]
tmp <- colSums(rce)/10^6
rcv2 <- t(t(rcv)/tmp)
if(sum(des[,"NormalizationGroup"]=="virus")==0){rcv2<-NULL}

rcb <- rawcnt[des[,"NormalizationGroup"]=="rRNA", ]
tmp <- colSums(rcb)/10^6
rcb2 <- t(t(rcb)/tmp)
if(sum(des[,"NormalizationGroup"]=="rRNA")==0){rcb2<-NULL}

rcc <- rawcnt[des[,"NormalizationGroup"]=="ercc", ]
tmp <- colSums(rcc)/10^6
rcc2 <- t(t(rcc)/tmp)
if(sum(des[,"NormalizationGroup"]=="ercc")==0){rcc2<-NULL}

rpm <- rbind(rce2, rcv2, rcb2, rcc2)
rpm[is.nan(rpm)] <- 0
rpm <- rpm[rownames(des),]


# save objects
save(rawcnt, file=sprintf("%s/%s", dir.input, fn.rawcnt))
save(rpm, file=sprintf("%s/%s", dir.input, fn.rpm))
save(cvrd1, file=sprintf("%s/%s", dir.input, fn.cvrd1))
save(cvrd3, file=sprintf("%s/%s", dir.input, fn.cvrd3))
save(readall, file=sprintf("%s/%s", dir.input, fn.readall))
save(readfw, file=sprintf("%s/%s", dir.input, fn.readfw))
save(fwpct, file=sprintf("%s/%s", dir.input, fn.fwpct))

write.table(cbind(readcount, at), file=sprintf("%s/%s", dir.output, fn.rcsummary))





