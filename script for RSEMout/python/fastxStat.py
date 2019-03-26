#!/usr/bin/env python
import sys, os
import gzip, bz2
from Bio import SeqIO
import numpy as np

headers = [
	"Sample name",
	"Total reads",
	"Total base",
	"Average length",
	"Median length",
	"Var of length",
	"Std of length",
	"Max length",
	"Min length"
]

def main():
	#========================================
	# parse args
	#----------------------------------------
	import argparse
	parser = argparse.ArgumentParser(description='This script outputs information of fasta(or fastq) file(s).')
	parser.add_argument("-i", metavar="fastq_file", dest="infiles", required=True, nargs="+", help='support gz file')
	parser.add_argument("-t", metavar="type", dest="type", required=True, choices=["fasta", "fastq"], help='fasta or fastq')
	args = parser.parse_args()
	
	print "\t".join(headers)
	
	all_lList = []
	#========================================
	# iterate files
	#----------------------------------------
	for infile in args.infiles:
		fileName = os.path.basename(infile)
		fileDir = os.path.dirname(infile)
		fileId, fileExt = os.path.splitext(fileName)
		
		#========================================
		# open file
		#----------------------------------------
		FHR = None
		if fileExt == ".gz":
			FHR = gzip.open(infile, "rb")
		elif fileExt == ".bz2":
			FHR = bz2.BZ2File(infile, "r")
		else:
			FHR = open(infile, "rU")
		
		#========================================
		# iterate reads
		#----------------------------------------
		lList = []
		for record in SeqIO.parse(FHR, args.type):
			id = record.id
			seq = record.seq
			l = len(seq)
			
			lList.append(l)
		FHR.close()
		
		all_lList.extend(lList)
		lList = np.array(lList)
		
		outputRow(fileId, lList)
	
	all_lList = np.array(all_lList)
	outputRow("Sum", all_lList)

def outputRow(rowname, lList):
	if len(lList) > 0:
		out = []
		out.append(rowname)
		out.append(len(lList))
		out.append(np.sum(lList))
		out.append(round(np.average(lList), 2))
		out.append(round(np.median(lList), 2))
		out.append(round(np.var(lList), 2))
		out.append(round(np.std(lList), 2))
		out.append(np.amax(lList))
		out.append(np.amin(lList))
		
		print "\t".join(map(str, out))
	else:
		print "\t".join(map(str, [rowname, 0, 0, 0, 0, 0, 0, 0, 0]))



# Run as script
if __name__=="__main__":
	main()


