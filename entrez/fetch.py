import os, sys
import pandas as pd
from Bio import SeqIO
from Bio import Entrez

if len(sys.argv) != 3:
	print("USAGE: python gb_fetch.py [ACCNO_LIST] [OUTPUT_FASTA]")
	sys.exit()

with open(sys.argv[1]) as handle:
    accnos = pd.Series([x.strip('" \n').lstrip('>') for x in handle.readlines()])

def split_groups(iterable, batch_size=5000):
	num_groups = len(accnos) // batch_size
	ranges = [(i*batch_size, (i+1)*batch_size) for i in range(num_groups)]
	ranges.append((num_groups * batch_size, len(iterable)))

	result = [iterable.iloc[start:end] for start, end in ranges]

	return result



print(f'Duplicated Count: {accnos.duplicated().sum()}')
accnos = accnos.drop_duplicates()
print(f'Final Count: {len(accnos)}')

Entrez.email = 'john.palmer1288@gmail.com'
Entrez.api_key = os.environ['NCBI_API_KEY']

if len(accnos) <= 5000:
	result = Entrez.efetch(db='nucleotide', id=accnos, rettype='fasta',retmode='text')
	with open(sys.argv[2], 'w') as outfile:
		for i in result:
			outfile.write(i)
else:
	with open(sys.argv[2], 'w') as outfile:
		for n, chunk in enumerate(split_groups(accnos)):
			print(f"Processing chunk {n}")
			result = Entrez.efetch(db='nucleotide', id=chunk.to_list(), rettype='fasta',retmode='text')
			for i in result:
				outfile.write(i)