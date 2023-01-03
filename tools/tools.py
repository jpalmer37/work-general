import os 
import sys
import pysam
from collections import defaultdict


#%%
def parse_fasta(filepath):	
	seqs = {}
	with open(filepath, 'r') as handle:
		for line in handle.readlines():
			if line[0] == '>':
				header = line.strip().lstrip('>')
				seqs[header] = ''
			else:
				seqs[header] += line.strip()
	return seqs
#%%
# helper function to handling read an input sam file
def alignment_reader(sam_input_path):
    """
    Function that returns a handle reading an input sam file
    """
    if sam_input_path:
        input_sam = pysam.AlignmentFile(sam_input_path, 'r')
    else:
        input_sam = pysam.AlignmentFile('-', 'r')
    return input_sam

# 
# def alignment_writer(bam_output_path, template):
#     """
#     Function that returns a handle that writes to a new output sam file
#     """
#     if bam_output_path:
#         output_bam = pysam.AlignmentFile(bam_output_path, 'wb', template=template)
#     else:
#         output_bam = pysam.AlignmentFile('-', 'wb', template=template)

#     return output_bam

def alignment_writer(bam_output_path, header):
    """
    Function that returns a handle that writes to a new output sam file
    """
    if bam_output_path:
        output_bam = pysam.AlignmentFile(bam_output_path, 'wb', header=header)
    else:
        output_bam = pysam.AlignmentFile('-', 'wb', header=header)

    return output_bam
#%%
def read_pair_generator(sam, region_string=None):
    """
    Generate read pairs for a SAM or BAM file or within a region string of said file.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in sam.fetch(region=region_string):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]

