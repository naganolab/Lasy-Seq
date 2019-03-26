################
### function ###
################

library("fields")
source("cb.R")
source("script/ng.Colors.R")
source("script/ng.BHFDR.R")


#set range
ng.mkrg <- function(input, mg.ratio=0.1){
  rg <- max(input)-min(input)
  mg <- rg*mg.ratio
  out <- c(min(input)-mg, max(input)+mg)
  return(out)
}


##############
### script ###
##############

# set parameters -----------------------------------------------------
dir.input <- "R"
dir.output <- "R"
dir.rawdata <- "result"

fn.description <- "141103-1_GeneDescription_Osa"
fn.attribute <- "SampleAttribute_Osa.txt"

fn.ercc <- "ERCC_Controls_Analysis.txt"

exec.date <- "190323"

fn.rawcnt <- sprintf("%s-2_rawcnt_Osa", exec.date)
fn.rpm <- sprintf("%s-3_rpm_Osa", exec.date)
fn.cvrd1 <- sprintf("%s-4_cvrd1_Osa", exec.date)
fn.cvrd3 <- sprintf("%s-5_cvrd3_Osa", exec.date)


#load sample attribute --------------
fn <- sprintf("%s/%s", dir.input, fn.attribute)
at <- read.delim(fn, header=T, as.is=T)
rownames(at) <- sprintf("%s_%03d", at[,"library"], at[,"index"])

#load Ahg transcript description --------------------
fn <- sprintf("%s/%s", dir.input, fn.description)
load(fn)

#set data read row --------------
drr <- des[,"NormalizationGroup"]=="data"
names(drr) <- rownames(des)

#load ERCC spike-in control annotation --------------
fn <- sprintf("%s/%s", dir.output, fn.ercc)
ercc.ano <- read.delim(fn, header=T, as.is=T)
rownames(ercc.ano) <- ercc.ano[,2]


# load expression data table-------------------------------------
load(file=sprintf("%s/%s", dir.input, fn.rawcnt))
load(file=sprintf("%s/%s", dir.input, fn.rpm))
load(file=sprintf("%s/%s", dir.input, fn.cvrd1))
load(file=sprintf("%s/%s", dir.input, fn.cvrd3))





# read count summay -------------------------------------------------
ng.ReadCountSummary <- function(rct, dr, ymax=100, m){
  hist(log10(colSums(rct)+1), breaks=seq(0,10,0.25), 
       xlab="log10(read count)", col="gray",
       main=sprintf("%s RNA-Seq", m))
  
  hist(log10(colSums(rct[dr,])+1), breaks=seq(0,10,0.25), 
       xlab="log10(read count)", col="gray",
       main=sprintf("%s RNA-Seq: data read", m))
  
  hist(log10(colSums(rct)+1), breaks=seq(0,10,0.25), 
       xlab="log10(read count)", ylim=c(0,ymax), col="gray",
       main=sprintf("%s RNA-Seq", m))
  
  hist(log10(colSums(rct[dr,])+1), breaks=seq(0,10,0.25), 
       xlab="log10(read count)", ylim=c(0,ymax), col="gray",
       main=sprintf("%s RNA-Seq: data read", m))
  
  plot(log10(colSums(rct)+1), log10(colSums(rct[dr,])+1), 
       xlim=c(0,9), ylim=c(0,9), xlab="mapped read", ylab="data read")
  segments(-1,-1,10,10, col="red")
  segments(0,-1,11,10, col="red")
  
  plot(log10(colSums(rct)), 
       main="mapped read")
  plot(log10(colSums(rct[dr,])), 
       main="data read")
  plot(colSums(rct[dr,])/colSums(rct),
       main="data read / mapped read")
  
}

pdf(file=sprintf("%s/190323-1_readcount.pdf", dir.output))

ng.ReadCountSummary(rawcnt, dr=drr, m="Osa readcount", ymax=100)

sample.set <- unique(at[,"set"])
for(k in sample.set){
  usecol <- (at[,"set"]==k)
  ng.ReadCountSummary(rawcnt[,usecol], dr=drr, 
                      m=sprintf("Arabidopsis halleri (%s)", k), 
                      ymax=100)
}

dev.off()


# make PerSampleHist -------------------------------------------------
pdf("R/190323-3_PerSampleHist.pdf")
tmp <- rawcnt[drr,]
rn <- colSums(tmp)
tmp <- tmp[,order(rn, decreasing=T)]
rn <- rn[order(rn, decreasing=T)]
tmp[log2(tmp)<(-2)] <- NA
for(i in 1:ncol(tmp)){
  hist(log2(tmp[, i]), 
       main=sprintf("%s, data read = %.0f", colnames(tmp)[i], rn[i]),
       breaks=seq(-2,23,1))
}
dev.off()


# make ERCC summary -------------------------------------------------
pdf("R/190323-5_ERCC summary.pdf")

tmp <- rawcnt[des[,"Type"]=="ercc_spikein",]
total.rn <- colSums(rawcnt)
ercc.conc <- ercc.ano[order(rownames(ercc.ano)),4]

plot(colSums(tmp)/total.rn, ylab="spikein / mapped read")
plot(colSums(tmp)/total.rn, ylab="spikein / mapped read", ylim=c(0,0.1))
plot(colSums(tmp)/total.rn, ylab="spikein / mapped read", ylim=c(0,0.01))
plot(sort(colSums(tmp)/total.rn), ylab="spikein / mapped read")
hist(colSums(tmp)/total.rn, xlab="spikein / mapped read", breaks=seq(0,1,0.025))

a <- cor(tmp, ercc.conc)
plot(a, ylab="r (rawcount)")
hist(a, xlab="r (rawcount)")
a <- cor(log2(tmp+1), log2(ercc.conc))
plot(a, ylab="r log2(rawcount+1)")
hist(a, xlab="r log2(rawcount+1)")

rn <- colSums(tmp)
tmp <- tmp[,order(rn, decreasing=T)]
rn <- rn[order(rn, decreasing=T)]
for(i in 1:ncol(tmp)){
  plot(log2(ercc.conc), log2(tmp[,i]+1), 
       main=sprintf("%s, ERCC spikein read = %.0f", colnames(tmp)[i], rn[i]))
}

dev.off()

tmp <- rawcnt[des[,"Type"]=="ercc_spikein",]
head(tmp)
rowSums(tmp)
a<-colSums(tmp)
head(a)
write.csv(a, "R/ERCCrawcont.csv")

pdf(file="R/190323ERCC contamination.pdf")
barplot(log10(rn), names.arg=names(rn),cex.names= 0.6,ylab="log10(rpm)",col=c("#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#888888","#888888","#888888","#888888"))
barplot(rn, names.arg=names(rn),cex.names= 0.6,ylab="log10(rpm)",col=c("#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#888888","#888888","#888888","#888888"))
barplot(rn[5:8], names.arg=names(rn[5:8]),cex.names= 0.6,ylab="log10(rpm)",col=c("#888888","#888888","#888888","#888888"))

dev.off()


# make rRNA summary -------------------------------------------------
pdf("R/190323-6_rRNAsummary.pdf")

tmp <- rawcnt[des[,"NormalizationGroup"]=="rRNA",]
total.rn <- colSums(rawcnt[des[,"NormalizationGroup"]!="ercc",])
plot(colSums(tmp)/total.rn, ylab="rRNA / mapped read")
plot(sort(colSums(tmp)/total.rn), ylab="rRNA / mapped read")

tmp <- rawcnt[des[,"NormalizationGroup"]=="rRNA",]
total.rn <- colSums(rawcnt)
ercc.conc <- ercc.ano[order(rownames(ercc.ano)),4]

plot(colSums(tmp)/total.rn, ylab="rRNA correlated / mapped read")
plot(sort(colSums(tmp)/total.rn), ylab="rRNA correlated / mapped read")
hist(colSums(tmp)/total.rn, xlab="rRNA correlated / mapped read", breaks=seq(0,1,0.025))

dev.off()


# make rpm+attribute+description table ----------------------------------------------

fn.RpmDesTable <- "190323_RpmDescription.csv"              
tmp <- rbind(t(at), rpm)
rpm.df <- data.frame(rpm)
rpm.df.des <- cbind(des, rpm.df)
fn <- sprintf("%s/%s", dir.output, fn.RpmDesTable)
write.csv(rpm.df.des, file=fn)

# make data read number table ----------------------------------------------

fn.DataReadNumber <- "190323_DataReadNumber.csv"
tmp <- rawcnt[drr, ]
sumDataReadNUm <- colSums(tmp)
fn <- sprintf("%s/%s", dir.output, fn.DataReadNumber)
write.csv(sumDataReadNUm, file=fn)




