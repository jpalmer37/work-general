import os, sys
import pandas as pd
from Bio import SeqIO
from Bio import Entrez

if len(sys.argv) != 3:
	print("USAGE: python gb_fetch.py [ACCNO_LIST] [OUTPUT_FASTA]")
	sys.exit()

with open(sys.argv[1]) as handle:
    accnos = pd.Series([x.strip('" \n').lstrip('>') for x in handle.readlines()])

print(f'Duplicated Count: {accnos.duplicated().sum()}')
accnos = accnos.drop_duplicates()
print(f'Final Count: {len(accnos)}')

Entrez.email = 'john.palmer1288@gmail.com'
Entrez.api_key = os.environ['NCBI_API_KEY']

result = Entrez.efetch(db='nucleotide', id=accnos, rettype='fasta',retmode='text')
with open(sys.argv[2], 'w') as outfile:
	for i in result:
		outfile.write(i)
