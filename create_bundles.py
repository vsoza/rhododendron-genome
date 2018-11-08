#this script takes a DAGchainer output with genomic coordinates from SynMap in CoGe and creates a bundle file for circos
#need to remove 1st comment line in DAGchainer output file before running python script

import sys
#open file at command line for reading
synfile = open(sys.argv[1], 'r')
#read file into a list of lines
synlist = synfile.readlines()
#close file
synfile.close()
#open new file for writing output via appending
bundlefile = open('bundles.txt', 'a')

#set counter
i=0
while i < len(synlist):
	#if line starts with #, extract number of syntenic genes in that syntenic block from end of line
	if synlist[i].startswith('#'):
		blocklist = synlist[i].split()
		pairs = int(blocklist[5])
		#then writes chromosome 1 and its start position from first pair, and end position from last pair
		#then writes chromoseom 2 and its start position from first pair, and end position from last pair
		pair1list = synlist[i+1].split()
		pair2list = synlist[i+pairs].split()
		bundlefile.write(pair1list[0]+'\t'+pair1list[2]+'\t'+pair2list[3]+'\t'+pair1list[4]+'\t'+pair1list[6]+'\t'+pair2list[7]+'\n')
		#then go to next line
		i=i+1
	else:
		i=i+1	
		
bundlefile.close()		
		