##############
### script ###
##############

# set parameters -----------------------------------------------------
dir.input <- "R"
dir.output <- "R"

fn.description <- "Araport11_genes201606transcript_rep_ERCC_Virus7457_desc"
fn.sample.index <- "190322_Sample-IndexPrimer.txt"
fn.sample.index.merged <- "190322_Sample-IndexPrimer_jmj_merged.txt"
fn.ercc <- "ERCC_Controls_Analysis.txt"

merge.sample <- TRUE
rpm.date <- "190322"
species <- "Ath"

fn.rawcnt <- sprintf("%s-2_rawcnt_%s", rpm.date, species)
fn.rpm <- sprintf("%s-3_rpm_%s", rpm.date, species)
fn.cvrd1 <- sprintf("%s-4_cvrd1_%s", rpm.date, species)
fn.cvrd3 <- sprintf("%s-5_cvrd3_%s", rpm.date, species)
fn.readall <- sprintf("%s-6_readall_%s", rpm.date, species)
fn.readfw <- sprintf("%s-7_readfw_%s", rpm.date, species)
fn.fwpct <- sprintf("%s-8_fwpct_%s", rpm.date, species)
fn.rcsummary <- sprintf("%s-1_ReadCountSummary_%s.txt", rpm.date, species)

exec.date <- "190322" 

fn.readcountsummary <- sprintf("%s-1_Summary_Readcount.pdf", exec.date)
fn.samplecorrelation <-  sprintf("%s-1_Summary_SampleCorrelation", exec.date)
fn.persamplehist <- sprintf("%s-1_Summary_PerSampleHist", exec.date)
fn.erccsummary <- sprintf("%s-1_Summary_ERCC", exec.date)
fn.rRNAsummary <- sprintf("%s-1_Summary_rRNA", exec.date)
fn.virussummary <- sprintf("%s-1_Summary_virus", exec.date)
fn.RpmDesTable <- sprintf("%s-1_RpmDescription", exec.date)
fn.DataReadNumber <- sprintf("%s-1_DataReadNumber", exec.date)


# load data ----------------------------------
#load sample attribute 
fn <- sprintf("%s/%s", dir.input, fn.sample.index)
if(merge.sample){fn <- sprintf("%s/%s", dir.input, fn.sample.index.merged)}
at <- read.delim(fn, header=T, as.is=T)
rownames(at) <- sprintf("%s_%03d", at[,"library"], at[,"index"])

#load transcript description 
fn <- sprintf("%s/%s", dir.input, fn.description)
load(fn)

#set data read row 
drr <- des[,"NormalizationGroup"]=="data"
names(drr) <- rownames(des)

#load ERCC spike-in control annotation 
fn <- sprintf("%s/%s", dir.input, fn.ercc)
ercc.ano <- read.delim(fn, header=T, as.is=T)
rownames(ercc.ano) <- ercc.ano[,2]

# load expression data table 
load(file=sprintf("%s/%s", dir.input, fn.rawcnt))
load(file=sprintf("%s/%s", dir.input, fn.rpm))
load(file=sprintf("%s/%s", dir.input, fn.cvrd1))
load(file=sprintf("%s/%s", dir.input, fn.cvrd3))
load(file=sprintf("%s/%s", dir.input, fn.readall))
load(file=sprintf("%s/%s", dir.input, fn.readfw))

# load expression read count summary table 
rcsummary <- read.table(file=sprintf("%s/%s", dir.input, fn.rcsummary))


# read count summay -------------------------------------------------
ng.ReadCountSummary <- function(rct, rcs, rall, rfwd, dr, ymax=NULL, m){
  
  tmp.hist <- function(...){
    hist(..., breaks=seq(0,10,0.25), 
         xlab="log10(read count)",
         ylim=ymax, 
         col="gray"
         )
  }
  
  tmp.hist(log10(rcs[,1]+1),
       main=sprintf("%s : total clean read", m))
  
  tmp.hist(log10(rcs[,2]+1),
       main=sprintf("%s : unique sequence", m))
  
  tmp.hist(log10(colSums(rct)+1),
       main=sprintf("%s : mapped read", m))
  
  tmp.hist(log10(colSums(rct[dr,])+1),
       main=sprintf("%s : data read", m))
  
  hist(100*colSums(rfwd[dr,])/colSums(rall[dr,]), breaks=seq(0,100,1),
       col="gray", main=sprintf("%s : strand specificity", m))
  
  plot(log10(rcs[,1]+1), log10(rcs[,2]+1), 
       xlim=c(0,9), ylim=c(0,9), xlab="total clean read", ylab="unique sequence")
  segments(-1,-1,10,10, col="red")
  segments(0,-1,11,10, col="red")
  
  plot(log10(colSums(rct)+1), log10(colSums(rct[dr,])+1), 
       xlim=c(0,9), ylim=c(0,9), xlab="mapped read", ylab="data read")
  segments(-1,-1,10,10, col="red")
  segments(0,-1,11,10, col="red")
  
  plot(log10(rcs[,1]+1), 
       main="total clean read")
  plot(log10(rcs[,2]+1), 
       main="unique sequence")
  plot(rcs[,2]/rcs[,1],
       main="unique sequence / total clean read")
  plot(log10(colSums(rct)), 
       main="mapped read")
  plot(colSums(rct)/rcs[,1],
       main="mapped read / total clean read")
  plot(log10(colSums(rct[dr,])), 
       main="data read")
  plot(colSums(rct[dr,])/colSums(rct),
       main="data read / mapped read")
  plot(colSums(rfwd[dr,])/colSums(rall[dr,]),
       main="strand specificity")
  
}

pdf(file=sprintf("%s/%s", dir.output, fn.readcountsummary))

ng.ReadCountSummary(rawcnt, rcs=rcsummary, rall=readall, rfwd=readfw, dr=drr, 
                    m="Read count summary")

sample.set <- unique(at[,"library"])
for(k in sample.set){
  usecol <- (at[,"library"]==k)
  ng.ReadCountSummary(rawcnt[,usecol], rcs=rcsummary[usecol,], 
                      rall=readall[,usecol], rfwd=readfw[,usecol], dr=drr, 
                      m=sprintf("Read count summary (%s)", k))
}
dev.off()


# make rpm+attribute+description table ----------------------------------------------

tmp <- rbind(t(at), rpm)

rpm.df <- data.frame(rpm)
rpm.df.des <- cbind(des, rpm.df)

fn <- sprintf("%s/%s.csv", dir.output, fn.RpmDesTable)
write.csv(rpm.df.des, file=fn)


# make data read number table ----------------------------------------------

tmp <- rawcnt[drr, ]
sumDataReadNUm <- colSums(tmp)

fn <- sprintf("%s/%s.csv", dir.output, fn.DataReadNumber)
DataReadNumber<-cbind(rownames(at),sumDataReadNUm)
write.csv(DataReadNumber, file=fn)

