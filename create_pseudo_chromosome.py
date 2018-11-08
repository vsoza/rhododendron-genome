#this script will take a multi-fasta file with a newline between fasta sequences and concatenate fasta sequences into one super-scaffold with 100 Ns between sequences

import sys
LG_file = open(sys.argv[1], 'r') 
#opens 'LG#_ordered.masked.fasta2' multi-fasta file as second command line argument
LG_list = LG_file.readlines()
#reads multi-fasta file into a list
LG_file.close()

pseudochromo_list = []
#creates new list for creating super-scaffold

#add each sequence from multi-fasta file to new list, without newline character, followed by 100 Ns, except for last sequence in multi-fasta file
i=1
while i < (len(LG_list)-1):
	pseudochromo_list.append(LG_list[i].rstrip()+("N" * 100))
	i = i + 3

pseudochromo_list.append(LG_list[-1])
#add last sequence from multi-fasta file to new list

sequence = "".join(pseudochromo_list)
#join sequence elements in list to create one super-sequence
pseudochromo_file = open('pseudochromo.fasta', 'w')
#create new fasta file for super-scaffold sequence
pseudochromo_file.write(LG_list[0])
#adds first fasta header from multi-fasta file to new fasta file
pseudochromo_file.write(sequence)
#adds super-sequence to new fasta file
pseudochromo_file.close()

