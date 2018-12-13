#script will extract GO terms for a csv list of genes from B2G .annot file

import sys

gene_file = open(sys.argv[1], 'r')
#opens csv file of gene names 

gene_string = gene_file.readline().strip()
#read 1st line of file as string and strip newline

gene_file.close()

gene_list = gene_string.split(',')
#turn gene string into a list

annot_file = open(sys.argv[2], 'r')
#opens B2G .annot file with gene names and GO terms

annot_list = []

for line in annot_file:
#for line in file, turn each line into a list, and append to annot_list
	line_list = line.split('\t')
	annot_list.append(line_list)

annot_file.close()

GO_file = open('GO_file.txt', 'a')
#create new file for appending

#for each gene in gene list, extract line from matching gene in .annot file and append to new GO file
for a in range(len(gene_list)):	
	for b in range(len(annot_list)):
		if gene_list[a] == annot_list[b][0]:
			gene_line = '\t'.join(annot_list[b])
			GO_file.write(gene_line)

GO_file.close()
