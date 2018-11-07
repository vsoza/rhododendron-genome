#modified code from rocky at https://stackoverflow.com/questions/31287065/how-to-filter-out-sequences-based-on-a-given-data-using-python
#filters out a set of sequences from a fasta file 

import sys
from Bio import SeqIO

remove_file = sys.argv[1] #fasta file of seqs to be deleted
fasta_file = sys.argv[2]  #original fasta file
output_file = sys.argv[3] #filtered fasta file

exclude = set()
fasta_sequences = SeqIO.parse(open(remove_file),'fasta')
for fasta in fasta_sequences:
    exclude.add(fasta.id)

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
with open(output_file, 'w') as output_handle:
   for fasta in fasta_sequences:
        if fasta.id not in exclude:
            SeqIO.write([fasta], output_handle, "fasta")