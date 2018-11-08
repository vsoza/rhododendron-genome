#For LG pseudochromosome, takes scaffold name and object end column from an agp file and creates an offset file for genome tools gt gff3 program that has scaffold name and offset value to use for adjusting feature positions in gff file
import sys

#open LG agp file 
#first line in file has scaffold name and running combined length of all previous scaffolds
#second line in file has running combined length of all previous scaffolds after adding 100 Ns to previous scaffold 
infile = open(sys.argv[1], 'r')

#read infile as a list
LG_agp_list = infile.readlines()
#open LG offset file for appending output from script to
outfile = open(sys.argv[2], 'a')

#write first scaffold to file with 0 offset
#split first line into a list at tab character
first_scaffold_list = LG_agp_list[0].split("\t")
#write scaffold name and length to outfile, formatted for genome tools offset file
outfile.write(first_scaffold_list[0]+"\t=\t0,\n") 

#set counters for scaffold name (i) and scaffold length (j)
i=2
j=1
while i < len(LG_agp_list) and j < len(LG_agp_list):
#split scaffold name line of file into list at tab character
	scaffold_list = LG_agp_list[i].split("\t")
#assign first list element to scaffold name
	scaffold_name = scaffold_list[0]
#write scaffold name to outfile, formatted for genome tools offset file
	outfile.write(scaffold_name+"\t=\t")
#assign running total from previous line to scaffold length, strip tab, newline characters, and return character from excel
	scaffold_length = LG_agp_list[j].strip()
#write scaffold length to outfile, formatted for genome tools offset file
	outfile.write(scaffold_length+",\n")
#go to next respective lines
	i=i+2
	j=j+2

#close files or else they wont be written to
infile.close()
outfile.close()	

#need to concatenate the resulting python files for all LGs and format end of file for genome tools
