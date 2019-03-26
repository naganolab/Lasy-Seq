#!/usr/bin/env python
import os, sys, re
import numpy as np


def main():
	import argparse
	parser = argparse.ArgumentParser(description='')
	parser.add_argument("-i", metavar="input.tsv", dest="input", required=True, nargs="+")
	parser.add_argument("-s", metavar="num of skip row", dest="skip", required=True, type=int)
	parser.add_argument("-c", metavar="id col index", dest="col", required=True, type=int)
	parser.add_argument("-v", metavar="value col index", dest="val", required=True, type=int)
	args = parser.parse_args()
	
	id_map = {}
	file_num = len(args.input)
	file_names = [""]
	
	file_index = 0
	for infile in args.input:
		fileName = os.path.basename(infile)
		fileDir = os.path.dirname(infile)
		fileId, fileExt = os.path.splitext(fileName)
		
		file_names.append(fileId)
		
		FHR = open(infile, "rU")
		for i in range(args.skip):
			FHR.next()
		for row in FHR:
			row = row.rstrip() # chomp
			itemList = row.split('\t')
			if len(itemList) < 1:
				continue
			
			id = itemList[args.col]
			if not id_map.has_key(id):
				id_map[id] = np.empty(file_num, dtype=object)
			id_map[id][file_index] = str(itemList[args.val])
			
		FHR.close()
		file_index += 1
	
	print "\t".join(file_names)
	
	for id in sorted(id_map.keys()):
		out = [id]
		for val in id_map[id]:
			if val == None:
				val = "."
			out.append(val)
		
		print "\t".join(out)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	# map_2 = {}
	
	# for infile in args.input2:
	# 	FHR = open(infile, "rU")
	# 	for i in range(args.skip2):
	# 		FHR.next()
	# 	for row in FHR:
	# 		row = row.rstrip() # chomp
	# 		itemList = row.split('\t')
	# 		if len(itemList) < 1:
	# 			continue
			
	# 		id = itemList.pop(args.col2)
	# 		if not map_2.has_key(id):
				
				
	# 			map_2[id] = itemList
	# 		# else:
	# 		# 	sys.exit("ERROR: ID duplication >> " + id)
	# 	FHR.close()
	
	
	
		
	# FHR = open(args.input1, "rU")
	# for i in range(args.skip1):
	# 	FHR.next()
	# for row in FHR:
	# 	row = row.rstrip() # chomp
	# 	itemList = row.split('\t')
	# 	if len(itemList) < 1:
	# 		continue
		
	# 	id = itemList[args.col1]
		
	# 	if map_2.has_key(id):
	# 		itemList.extend(map_2[id])
		
	# 	print "\t".join(itemList)
		
	# FHR.close()
	
# Run as script
if __name__=="__main__":
	main()




