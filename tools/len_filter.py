import os, sys
from Bio import SeqIO


if len(sys.argv) != 4:
	print("USAGE: python len_filter.py [INPUT_FASTA] [MIN_LENGTH] [OUTPUT_FASTA]")
	sys.exit()

seqs = SeqIO.to_dict(SeqIO.parse(sys.argv[1],'fasta'))

filtered = {x:y for x,y in seqs.items() if len(y) > int(sys.argv[2])}

SeqIO.write(list(filtered.values()), sys.argv[3], 'fasta')