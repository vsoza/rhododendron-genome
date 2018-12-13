#script will extract functions from a Blast2GO annotation file for syntenic genes within a genome identified by SynMap2
#you will need the B2G file in annot format, gff file, and modified DAGchainer output file, see below

import sys

DAGchainer_file = open(sys.argv[1], 'r')
#opens modified DAGchainer output file that has pseudochromosome name and genomic coordinates of queries

DAGchainer_list=[]
#create new DAGchainer list

for line in DAGchainer_file:
#for line in file, turn each line into a list, and append to DAGchainer list
	line_list = line.split()
	DAGchainer_list.append(line_list)

DAGchainer_file.close()

gff_file = open(sys.argv[2], 'r')
#opens gff file for pseudochromosomes

gff_list=[]
#create new gff_list

for line in gff_file:
#for line in file, turn each line into a list, and append to gff list
	line_list = line.split()
	gff_list.append(line_list)

gff_file.close()

gff_col9_list = []
#create new list for adding column 9 from gff file for matches to DAGchainer file

#for each query in DAGchainer list, extract column 9 from matching mRNA in gff file and append to gff_col9_list
#if name of pseudochromosome is the same in DAGchainer file and gff file, and feature is an mRNA, and start and stop positions are the same (including reverse orientation), extract column 9
for h in range(len(DAGchainer_list)):	
	for i in range(len(gff_list)):
		if (DAGchainer_list[h][0] == gff_list[i][0]) and (gff_list[i][2] == "mRNA") and (DAGchainer_list[h][1] == gff_list[i][3] or DAGchainer_list[h][1] == gff_list[i][4]) and (DAGchainer_list[h][2] == gff_list[i][3] or DAGchainer_list[h][2] == gff_list[i][4]):
			gff_col9_list.append(gff_list[i][8])

gff_col9_list2 = []
#create new list for appending each (gff) column 9 as a list, see below

#turn each column 9 into a list
for item in gff_col9_list:
	item_list = item.split(";")
	gff_col9_list2.append(item_list)

name_list = []
#create new list for appending only the name field from each (gff) column 9, see below	

#extract name field from each column 9 and add to name_list
j=0
while j < len(gff_col9_list2):
	name_list.append(gff_col9_list2[j][2])
	j=j+1
	
name_list2 = []
#create new list for appending only the name value from below

#turn name field into a list and add name value to name_list2
for item in name_list:
	item_list = item.split("=")
	name_list2.append(item_list[1].strip())

syngene_annot_file = open('syntenic_gene.annot', 'a')
#open file for appending

B2G_file = open(sys.argv[3], 'r')
#opens Blast2GO annotation file exported as annot format, which has name of protein and 1 GO term plus sequence description per line
B2G_list = B2G_file.readlines()
#reads B2G file as a list of lines
B2G_file.close()

B2G_list2 = []
#create new list for B2G annotations

#turn B2G annotations lines into lists
for item in B2G_list:
	item_list = item.split("\t")
	B2G_list2.append(item_list)

#for each gene in name_list2, extract annotation from B2G file and write to new file
for l in range(len(name_list2)):
	for k in range(len(B2G_list2)):
		if name_list2[l] == B2G_list2[k][0]:
			syngene_annot_file.write(B2G_list[k])
			 
syngene_annot_file.close()
