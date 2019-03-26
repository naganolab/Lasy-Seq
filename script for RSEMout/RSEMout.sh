#!/usr/bin/env bash
set -eu
BINDIR=/opt/bio/work/bin

## working directory
BASEDIR=/home/kamitani/jmjKM2_180502/1804KHP-0020
## reference
DB=/home/kamitani/jmjKM2_180502/1804KHP-0020/data/Araport11/Araport11_genes.201606.transcript.rep_ERCC_Virus7457.fa
## bed file
BED=${DB}.bed
#BED=/home/kamitani/jmjKM2_180502/1804KHP-0020/data/Araport11/Araport11_genes.201606.cdna.rep_ERCC_Virus7457.bed  


echo "### makeRSEMout ========================================================="
$BINDIR/batch_makeRSEMout_noid.sh \
	$OUTDIR \
	$FINALDIR \
	$BED
date

if [ ! -e $BASEDIR/data/count_reads.tsv ]; then
	FILES=`ls $BASEDIR/data/org/*.fastq.gz; ls $BASEDIR/data/preprocessed/*.fq`
	echo "### count reads"
	/opt/bio/work/python/fastxStat.py -t fastq -i $FILES > $BASEDIR/data/count_reads.tsv
	date
fi

if [ ! -e $FINALDIR/count_mapped_reads.tsv ]; then
	echo "### count mapped reads"
	FILES=`ls $FINALDIR/*.results`
	/opt/bio/work/python/parseRSEMresult.py -i $FILES > $FINALDIR/count_mapped_reads.tsv
	date
fi

if [ ! -e $FINALDIR/merged_stat.tsv ]; then
	echo "### merge TSVs"
	cat $FINALDIR/count_mapped_reads.tsv |grep -v "^#"|  cut -f 1,4,5 | sed -e "s/.genes//" > $FINALDIR/count_mapped_reads.shaped.tsv
	/opt/bio/work/python/mergeById.py -1 $BASEDIR/data/count_reads.tsv -c1 0 -s1 1 -2 $FINALDIR/count_mapped_reads.shaped.tsv -c2 0 -s2 0 > $FINALDIR/merged_stat.tsv
	date
fi


rm -f $OUTDIR/*.transcript.bam

# items=$(ls $BASEDIR/data/preprocessed/*.fq)
# for item in ${items[@]}; do
# 	gzip $item
# done
