#!/usr/bin/env python
import sys, os, re
from subprocess import Popen, PIPE, check_call
import numpy as np

SAMTOOLS = "samtools"

def main():
	#========================================
	# parse args
	#----------------------------------------
	import argparse
	parser = argparse.ArgumentParser(description='This script outputs information of bam file(s).')
	parser.add_argument("-i", metavar="RSEM_result_file", dest="infiles", required=True, nargs="+")
	args = parser.parse_args()
	
	#========================================
	# iterate files
	#----------------------------------------
	c = 0
	for infile in args.infiles:
		fileName = os.path.basename(infile)
		fileDir = os.path.dirname(infile)
		fileId, fileExt = os.path.splitext(fileName)
		
		#========================================
		# open file
		#----------------------------------------
		cmd = "cat " + infile + " | grep '^#'"
		p = Popen([cmd], shell=True, stdout=PIPE)
		
		#========================================
		# iterate records
		#----------------------------------------
		titles = ["sample name"]
		values = [fileId]
		for row in p.stdout:
			row = row.rstrip() # chomp
			itemList = row.split('\t')
			if len(itemList) < 1:
				continue
			
			title = re.sub("^#", "", itemList[0])
			titles.append(title)
			values.append(itemList[1])
		p.wait()
		p.stdout.close()
		
		if c < 1:
			print "\t".join(titles)
		print "\t".join(values)
		
		
		c += 1


# Run as script
if __name__=="__main__":
	main()


