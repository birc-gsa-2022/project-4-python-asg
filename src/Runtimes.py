################################################################
# libraries:
import re
import time
import random
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

################################################################
# functions:
    
from SEQsimulator import simulate_string
from SEQsimulator import get_exact_read


def re_find(string, pattern):
    re_findings = [m.start() for m in re.finditer('(?={0})'.format(re.escape(pattern)), string)]
    return re_findings 


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


def get_features(string):
    string+='$'
    alphabet = sorted(set(string))
    counts = {v: string.count(v) for v in alphabet}
    index = {v: i for i, v in enumerate(alphabet)}
    return alphabet, counts, index


def burrow_wheeler(string):
    SA = SuffixArray(string)
    F = [None]*len(SA)
    L = F[:]
    reduced_SA = {}
    string += '$'
    for i in range(len(SA)):
        F[i] = string[SA[i]]
        L[i] = string[SA[i]-1]
        if SA[i]%10 == 0:
            reduced_SA[i] = SA[i]
    return reduced_SA, F, L


def reduced_L_tally(L, alphabet):
    counts = {v: 0 for v in alphabet}
    tally = {}
    for i in range(len(L)):
        counts[L[i]]+=1
        if i%5 == 0 or i==len(L)-1:
            tally[i] = tuple(counts.values())
    return tally


def FM_structures(string):
    alphabet, quant, index = get_features(string)
    red_SA, F, L = burrow_wheeler(string)
    red_L_tally = reduced_L_tally(L, alphabet)
    return alphabet, quant, index, red_SA, F, L, red_L_tally


def get_L_interval(letter: str):
    start = 0
    for i in range(len(alphabet)-1):
        if alphabet[i] == letter:
            break
        start += quant[alphabet[i]]
    end = start + quant[letter]
    interval = (start, end)
    return interval


def get_L_rows(letter: str, interval: tuple):
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


def LF_mapping(letter, n, first_rank):
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


def find_pattern_interval(pattern):
    j = len(pattern)-1
    interval = get_L_interval(pattern[j])
    while j > 0:
        n, first_rank = get_L_rows(pattern[j-1], interval)
        interval = LF_mapping(pattern[j-1], n, first_rank)
        j-=1
    if interval[0] == interval[1]:
        return []
    return interval


def get_SA_offsets(SA_interval: tuple):
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
                L_row = LF_mapping(letter, 1, rank)[0]
                if L_row in red_SA:
                    SA_val = red_SA[L_row] + steps
                    offsets.append(SA_val)
                    S=1
                else:
                    k = L_row
    return offsets

################################################################
# tests:
  
# 1 mio test:
string = simulate_string(1000000)
start_time = time.time()
FM_structures(string)
end_time = time.time()
print(end_time-start_time)


# Test suffixtree-algorithm vs naive-algorithm for same result:
for i in range(50000):
    print('Iteration nr: ', i+1)
    ref = simulate_string(random.randint(40,100))
    read = get_exact_read(ref, random.randint(1,30))
    alphabet, quant, index, red_SA, F, L, red_L_tally = FM_structures(ref)
    SA_interval = find_pattern_interval(read)
    offsets = get_SA_offsets(SA_interval)
    if sorted(offsets) != sorted(re_find(ref,read)):
        print('Algorithm mistake!')
        print(ref)
        print(read)
        break
    if i+1 == 50000:
        print('DONE')


# Runtimes for the array construction (varying ref lengths):
ref_lengths = [25000,50000,75000,100000,125000,150000,175000]
runtimes = []
for idx in range(7):
    print('Iteration nr: ', idx+1) 
    replicate = []
    for j in range(10):
        ref = simulate_string(ref_lengths[idx])
        start_time = time.time()
        FM_structures(ref)
        end_time = time.time()
        replicate.append(end_time-start_time)
    runtimes.append(np.mean(replicate))
# plot running times:
fig, ax = plt.subplots()
sns.lineplot(x=ref_lengths, y=runtimes, ax=ax)
plt.xlabel('ref length')
plt.ylabel('runtime (s)')
plt.tight_layout()
plt.show()


# Runtimes for mapping (varying read lengths):
ref_lengths = [2500,5000,7500,10000,12500]
read_lengths_10 = [100]*5
read_lengths_20 = [200]*5
read_lengths_30 = [300]*5
read_lengths_40 = [400]*5
read_lengths_50 = [500]*5
runtimes_10 = []
runtimes_20 = []
runtimes_30 = []
runtimes_40 = []
runtimes_50 = []
for idx in range(5):
    print('Iteration nr: ', idx+1)
    runtimes_10_replicate = []
    runtimes_20_replicate = []
    runtimes_30_replicate = []
    runtimes_40_replicate = []
    runtimes_50_replicate = []
    
    for i in range(10):
        ref = simulate_string(ref_lengths[idx])
        alphabet, quant, index, red_SA, F, L, red_L_tally = FM_structures(ref)
        
        read_10 = get_exact_read(ref, read_lengths_10[idx])
        read_20 = get_exact_read(ref, read_lengths_20[idx])
        read_30 = get_exact_read(ref, read_lengths_30[idx])
        read_40 = get_exact_read(ref, read_lengths_40[idx])
        read_50 = get_exact_read(ref, read_lengths_50[idx])
        
        # Dont know why it is nessesary to read run this chunk of code, but 
        # if i dont the first the first runtime is affected. Mayde something 
        # to do with loading the modules??
        SA_interval = find_pattern_interval(read_10)
        offsets = get_SA_offsets(SA_interval)
        
        start_time = time.time()
        SA_interval = find_pattern_interval(read_10)
        offsets = get_SA_offsets(SA_interval)
        end_time = time.time()
        runtimes_10_replicate.append(end_time-start_time)
        
        start_time = time.time()
        SA_interval = find_pattern_interval(read_20)
        offsets = get_SA_offsets(SA_interval)
        end_time = time.time()
        runtimes_20_replicate.append(end_time-start_time)
        
        start_time = time.time()
        SA_interval = find_pattern_interval(read_30)
        offsets = get_SA_offsets(SA_interval)
        end_time = time.time()
        runtimes_30_replicate.append(end_time-start_time)
        
        start_time = time.time()
        SA_interval = find_pattern_interval(read_40)
        offsets = get_SA_offsets(SA_interval)
        end_time = time.time()
        runtimes_40_replicate.append(end_time-start_time)
        
        start_time = time.time()
        SA_interval = find_pattern_interval(read_50)
        offsets = get_SA_offsets(SA_interval)
        end_time = time.time()
        runtimes_50_replicate.append(end_time-start_time)
        
    runtimes_10.append(np.mean(runtimes_10_replicate))
    runtimes_20.append(np.mean(runtimes_20_replicate))
    runtimes_30.append(np.mean(runtimes_30_replicate))
    runtimes_40.append(np.mean(runtimes_40_replicate))
    runtimes_50.append(np.mean(runtimes_50_replicate))

# plot running times:
fig, ax = plt.subplots()
sns.lineplot(x=ref_lengths, y=runtimes_50, ax=ax, label='read length = 50')
sns.lineplot(x=ref_lengths, y=runtimes_40, ax=ax, label='read length = 40')
sns.lineplot(x=ref_lengths, y=runtimes_30, ax=ax, label='read length = 30')
sns.lineplot(x=ref_lengths, y=runtimes_20, ax=ax, label='read length = 20')
sns.lineplot(x=ref_lengths, y=runtimes_10, ax=ax, label='read length = 10')
plt.xlabel('ref length')
plt.ylabel('runtime (s)')
plt.tight_layout()
plt.show()

################################################################
