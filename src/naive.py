"""Implementation of the naive exact matching algorithm."""

################################################################
# libraries:
import sys

################################################################
# Functions:

from align import get_edits
from cigar import edits_to_cigar    

def naive_algorithm(ref:str, read:str):
    if ref == '' or ref == None:
        return []
    if read == '' or read == None:
        return []
    ref=ref.strip().upper()
    read=read.strip().upper()
    if len(read) > len(ref):
        return False
    if len(read) == len(ref):
        return False
    
    ref=ref.strip().upper()
    read=read.strip().upper()
    match_positions = []
    for idx in range(len(ref)-len(read)+1):
        substring = ref[idx:idx+len(read)]
        for i in range(len(substring)):
            if substring[i] != read[i]:
                break
            if i == len(substring)-1:
                if substring[i] == read[i]:
                    match_positions.append(idx)
    return match_positions


def read_fasta():
    # load input:
    inFile = sys.argv[1]
    with open(inFile,'r') as f:
        lines = f.readlines()
    record_list = []
    header = ''
    sequence = []
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            if header != "":
                record_list.append([header.strip(), ''.join(sequence).strip()])
                sequence = []
            header = line[1:]
        else:
            sequence.append(line)
    record_list.append([header.strip(), ''.join(sequence).strip()])
    return record_list


def read_fastq():
    inFile = sys.argv[2]  
    with open(inFile,'r') as f:
        lines = f.readlines()
    record_list = []
    header = ''
    sequence = []
    for line in lines:
        line = line.strip()
        if line.startswith('@'):
            if header != "":
                record_list.append([header.strip(), ''.join(sequence).strip()])
                sequence = []
            header = line[1:]
        else:
            sequence.append(line)
    record_list.append([header.strip(), ''.join(sequence).strip()])
    return record_list

################################################################
# Code:
    
if __name__ == '__main__':
    
    fasta_recs = read_fasta()
    fastq_recs = read_fastq()
    
    for fq_rec in fastq_recs:
        for fa_rec in fasta_recs:
            matches = naive_algorithm(fa_rec[1], fq_rec[1])
            for match in matches:
                read_name = fq_rec[0]
                read_seq = fq_rec[1]
                edits = get_edits(read_seq, fa_rec[1][match:match+len(fq_rec[1])])
                cigar = edits_to_cigar(edits[2])
                output = [read_name,fa_rec[0],str(match+1),cigar,read_seq]
                print('\t'.join(output))
        

################################################################