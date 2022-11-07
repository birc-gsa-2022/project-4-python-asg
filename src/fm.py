####################################################
# Usage: 
# string = 'mississippi'
# pattern = ''
# alphabet, quant, index, red_SA, F, L, red_L_tally = FM_structures(string)
# SA_interval = find_pattern_interval(pattern, alphabet, index, red_L_tally, quant, L)
# offsets = get_SA_offsets(SA_interval, red_SA, L, red_L_tally, index, alphabet, quant)
# print(offsets)

####################################################
# Libraries: 
    
import sys
import argparse
import os

####################################################
# Functions: 
    
from cigar import edits_to_cigar
from align import get_edits
    
def SuffixArray(string):
    if string == '' or string == None:
        return None
    string += '$'
    index = {v: i for i, v in enumerate(sorted(set(string)))}
    string = [index[v] for v in string]
    rank_list = list(range(len(string)))
    SA = [None]*len(string)
    tuple_list = SA[:]
    M,j = 0,1
    while M < len(string)-1:
        for i in range(len(string)):
            i_j = i+j
            if i_j < len(string):
                key = (string[i], string[i_j])
            else: 
                key = (string[i], string[-1])
            tuple_list[i] = key
        j*=2
        keys = sorted(set(tuple_list))
        ranks = dict(zip(keys, rank_list))
        string = [ranks[tuple_list[i]] for i in range(len(string))]
        M = max(string)
    for i in rank_list:
        SA[string[i]] = i
    return SA


def FM_structures(string):
    if string == '' or string == None or string == []:
        return []
    
    SA = SuffixArray(string)
    
    # Get features:
    string+='$'
    alphabet = sorted(set(string))
    quant = {v: string.count(v) for v in alphabet}
    index = {v: i for i, v in enumerate(alphabet)}
    
    # BWT:
    F = [None]*len(SA)
    L = F[:]
    reduced_SA = {}
    string += '$'
    for i in range(len(SA)):
        F[i] = string[SA[i]]
        L[i] = string[SA[i]-1]
        if SA[i]%10 == 0:
            reduced_SA[i] = SA[i]
    
    # Make reduced table:
    counts = {v: 0 for v in alphabet}
    tally = {}
    for i in range(len(L)):
        counts[L[i]]+=1
        if i%5 == 0 or i==len(L)-1:
            tally[i] = tuple(counts.values())
    
    red_SA = reduced_SA
    red_L_tally = tally
    return alphabet, quant, index, red_SA, F, L, red_L_tally


def get_L_interval(letter: str, alphabet, quant):
    start = 0
    for i in range(len(alphabet)-1):
        if alphabet[i] == letter:
            break
        start += quant[alphabet[i]]
    end = start + quant[letter]
    interval = (start, end)
    return interval


def get_L_rows(letter, interval, index, red_L_tally, L):
    above_val = 0 
    end_val = 0 
    counts = 0
    for i in range(interval[0]-1, -1, -1):
        if i in red_L_tally:
            above_val = red_L_tally[i][index[letter]] + counts
            break
        if L[i] == letter:
            counts+=1
    counts = 0
    for i in range(interval[1], len(L)):
        if L[i] == letter:
            counts+=1
        if i in red_L_tally:
            end_val = red_L_tally[i][index[letter]] - counts
            break
    if interval[1] == len(L): 
        end_val = red_L_tally[len(L)-1][index[letter]]
    n = end_val - above_val
    rank_first = above_val
    return n, rank_first


def LF_mapping(letter, n, first_rank, alphabet, quant):
    accumulator = 0
    for i in range(len(alphabet)):
        if alphabet[i] != letter:
            accumulator += quant[alphabet[i]]
        else:
            accumulator += first_rank
            break
    start = accumulator
    end = accumulator + n
    return (start, end)


def find_pattern_interval(pattern, alphabet, index, red_L_tally, quant, L):
    if pattern == '' or pattern == None or pattern == []:
        return []
        
    j = len(pattern)-1
    interval = get_L_interval(pattern[j], alphabet, quant)
    while j > 0:
        n, first_rank = get_L_rows(pattern[j-1], interval, index, red_L_tally, L)
        interval = LF_mapping(pattern[j-1], n, first_rank, alphabet, quant)
        j-=1
    if interval[0] == interval[1]:
        return []
    return interval


def get_SA_offsets(SA_interval, red_SA, L, red_L_tally, index, alphabet, quant):
    
    if SA_interval == '' or SA_interval == None or SA_interval == []:
        return []
    
    offsets = []
    for i in range(SA_interval[0], SA_interval[1]):
        steps = 0
        if i in red_SA: 
            offsets.append(red_SA[i])
            continue
        else:
            S=0
            k=i
            while S==0:
                steps += 1
                # get L-rank
                rank = None
                counts = 0
                letter = L[k]
                for j in range(k-1, -1, -1): 
                    if j in red_L_tally:
                        rank = red_L_tally[j][index[letter]] + counts
                        break
                    if L[j] == letter:
                        counts+=1   
                # get new L-row
                L_row = LF_mapping(letter, 1, rank, alphabet, quant)[0]
                if L_row in red_SA:
                    SA_val = red_SA[L_row] + steps
                    offsets.append(SA_val)
                    S=1
                else:
                    k = L_row
    return offsets


def read_fasta(inFile):
    lines = inFile.readlines()
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


def read_fastq(inFile):
    lines = inFile.readlines()
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


def write_FM_structures(genome_name, name, alphabet, quant, index, red_SA, F, L, red_L_tally):
    
    os.makedirs('{}/{}/'.format(genome_name,name))
    
    with open('{}/{}/alphabet.txt'.format(genome_name,name), 'w') as f:
        print(alphabet, file=f)
        
    with open('{}/{}/quant.txt'.format(genome_name,name), 'w') as f:
        print(quant, file=f)
    
    with open('{}/{}/index.txt'.format(genome_name,name), 'w') as f:
        print(index, file=f)
    
    with open('{}/{}/red_SA.txt'.format(genome_name,name), 'w') as f:
        print(red_SA, file=f)
        
    with open('{}/{}/F.txt'.format(genome_name,name), 'w') as f:
        print(F, file=f)
    
    with open('{}/{}/L.txt'.format(genome_name,name), 'w') as f:
        print(L, file=f)
        
    with open('{}/{}/red_L_tally.txt'.format(genome_name,name), 'w') as f:
        print(red_L_tally, file=f)
        

def open_FM_structures(path_to_dir: str):
    alphabet = open('alphabet.txt', 'r').read()
    alphabet = eval(alphabet)
    
    quant = open('quant.txt', 'r').read()
    quant = eval(quant)
    
    index = open('index.txt', 'r').read()
    index = eval(index)
    
    red_SA = open('red_SA.txt', 'r').read()
    red_SA = eval(red_SA)
    
    F = open('F.txt', 'r').read()
    F = eval(F)
    
    L = open('L.txt', 'r').read()
    L = eval(L)
    
    red_L_tally = open('red_L_tally.txt', 'r').read()
    red_L_tally = eval(red_L_tally)
    
    return alphabet, quant, index, red_SA, F, L, red_L_tally

####################################################
# Main:

def main():
    argparser = argparse.ArgumentParser(
        description="FM-index exact pattern matching",
        usage="\n\tfm -p genome\n\tfm genome reads"
    )
    argparser.add_argument(
        "-p", action="store_true",
        help="preprocess the genome."
    )
    argparser.add_argument(
        "genome",
        help="Simple-FASTA file containing the genome.",
        type=argparse.FileType('r')
    )
    argparser.add_argument(
        "reads", nargs="?",
        help="Simple-FASTQ file containing the reads.",
        type=argparse.FileType('r')
    )
    args = argparser.parse_args()

    if args.p:
        # print(f"Preprocess {args.genome}")
        fasta_recs = read_fasta(args.genome)
        for fa_rec in fasta_recs:
            alphabet, quant, index, red_SA, F, L, red_L_tally = FM_structures(fa_rec)
            try: genome_name = args.genome.name.split('/')[1]
            except: genome_name = args.genome.name
            write_FM_structures('{}'.format(genome_name),'Preprocessed_{}'.format(fa_rec[0]), alphabet, quant, index, red_SA, F, L, red_L_tally)
    else:
        # print(f"Search {args.genome} for {args.reads}")
        if args.reads is None:
            argparser.print_help()
            sys.exit(1)
        else:
            fasta_recs = read_fasta(args.genome)
            fastq_recs = read_fastq(args.reads)
            print(fasta_recs)
            for fa_rec in fasta_recs:
                ref = fa_rec[1]
                alphabet, quant, index, red_SA, F, L, red_L_tally = FM_structures(ref)
                for fq_rec in fastq_recs:
                    read = fq_rec[1]
                    SA_interval = find_pattern_interval(read, alphabet, index, red_L_tally, quant, L)
                    matches = get_SA_offsets(SA_interval, red_SA, L, red_L_tally, index, alphabet, quant)
                    for match in matches:
                        read_name = fq_rec[0]
                        read_seq = fq_rec[1]
                        edits = get_edits(read_seq, fa_rec[1][match:match+len(fq_rec[1])])
                        cigar = edits_to_cigar(edits[2])
                        output = [read_name,fa_rec[0],str(match+1),cigar,read_seq]
                        print('\t'.join(output))
        
####################################################
# Code: 

if __name__ == '__main__':
    main()
    
####################################################
