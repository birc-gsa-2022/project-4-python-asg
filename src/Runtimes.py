################################################################
# libraries:
import re
import time
import random
# import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns

################################################################
# functions:
    
from SEQsimulator import simulate_string
from SEQsimulator import get_exact_read
from fm import FM_structures
from fm import find_pattern_interval
from fm import get_SA_offsets
from naive import naive_algorithm


def re_find(string, pattern):
    re_findings = [m.start() for m in re.finditer('(?={0})'.format(re.escape(pattern)), string)]
    return re_findings 


def FM_read_mapper(string, pattern):
    alphabet, quant, index, red_SA, F, L, red_L_tally = FM_structures(string)
    SA_interval = find_pattern_interval(pattern, alphabet, index, red_L_tally, quant, L)
    offsets = get_SA_offsets(SA_interval, red_SA, L, red_L_tally, index, alphabet, quant)
    return offsets

################################################################
# # tests:
  
# 1 mio test:
# string = simulate_string(1000000)
# start_time = time.time()
# FM_structures(string)
# end_time = time.time()
# print(end_time-start_time)


# Test suffixtree-algorithm vs naive-algorithm for same result:
# for i in range(50000):
#     print('Iteration nr: ', i+1)
#     ref = simulate_string(random.randint(40,100))
#     read = get_exact_read(ref, random.randint(0,30))
#     offsets = FM_read_mapper(ref, read)
#     if sorted(offsets) != sorted(naive_algorithm(ref,read)):
#         print('Algorithm mistake!')
#         print(ref)
#         print(read)
#         break
#     if i+1 == 50000:
#         print('DONE')


# Runtimes for the array construction (varying ref lengths):
# ref_lengths = [25000,50000,75000,100000,125000,150000,175000]
# runtimes = []
# for idx in range(7):
#     print('Iteration nr: ', idx+1) 
#     replicate = []
#     for j in range(10):
#         ref = simulate_string(ref_lengths[idx])
#         start_time = time.time()
#         FM_structures(ref)
#         end_time = time.time()
#         replicate.append(end_time-start_time)
#     runtimes.append(np.mean(replicate))
# # plot running times:
# fig, ax = plt.subplots()
# sns.lineplot(x=ref_lengths, y=runtimes, ax=ax)
# plt.xlabel('ref length')
# plt.ylabel('runtime (s)')
# plt.tight_layout()
# plt.show()


# Runtimes for mapping (varying read lengths):
# ref_lengths = [2500,5000,7500,10000,12500]
# read_lengths_10 = [10]*5
# read_lengths_20 = [20]*5
# read_lengths_30 = [30]*5
# read_lengths_40 = [40]*5
# read_lengths_50 = [50]*5
# runtimes_10 = []
# runtimes_20 = []
# runtimes_30 = []
# runtimes_40 = []
# runtimes_50 = []
# for idx in range(5):
#     print('Iteration nr: ', idx+1)
#     runtimes_10_replicate = []
#     runtimes_20_replicate = []
#     runtimes_30_replicate = []
#     runtimes_40_replicate = []
#     runtimes_50_replicate = []
    
#     for i in range(100):
#         ref = simulate_string(ref_lengths[idx])
#         alphabet, quant, index, red_SA, F, L, red_L_tally = FM_structures(ref)
        
#         read_10 = get_exact_read(ref, read_lengths_10[idx])
#         read_20 = get_exact_read(ref, read_lengths_20[idx])
#         read_30 = get_exact_read(ref, read_lengths_30[idx])
#         read_40 = get_exact_read(ref, read_lengths_40[idx])
#         read_50 = get_exact_read(ref, read_lengths_50[idx])
        
#         # Dont know why it is nessesary to read run this chunk of code, but 
#         # if i dont the first the first runtime is affected. Mayde something 
#         # to do with loading the modules??
#         SA_interval = find_pattern_interval(read_10, alphabet, index, red_L_tally, quant, L)
#         offsets = get_SA_offsets(SA_interval, red_SA, L, red_L_tally, index, alphabet, quant)
        
#         start_time = time.time()
#         SA_interval = find_pattern_interval(read_10, alphabet, index, red_L_tally, quant, L)
#         offsets = get_SA_offsets(SA_interval, red_SA, L, red_L_tally, index, alphabet, quant)
#         end_time = time.time()
#         runtimes_10_replicate.append(end_time-start_time)
        
#         start_time = time.time()
#         SA_interval = find_pattern_interval(read_20, alphabet, index, red_L_tally, quant, L)
#         offsets = get_SA_offsets(SA_interval, red_SA, L, red_L_tally, index, alphabet, quant)
#         end_time = time.time()
#         runtimes_20_replicate.append(end_time-start_time)
        
#         start_time = time.time()
#         SA_interval = find_pattern_interval(read_30, alphabet, index, red_L_tally, quant, L)
#         offsets = get_SA_offsets(SA_interval, red_SA, L, red_L_tally, index, alphabet, quant)
#         end_time = time.time()
#         runtimes_30_replicate.append(end_time-start_time)
        
#         start_time = time.time()
#         SA_interval = find_pattern_interval(read_40, alphabet, index, red_L_tally, quant, L)
#         offsets = get_SA_offsets(SA_interval, red_SA, L, red_L_tally, index, alphabet, quant)
#         end_time = time.time()
#         runtimes_40_replicate.append(end_time-start_time)
        
#         start_time = time.time()
#         SA_interval = find_pattern_interval(read_50, alphabet, index, red_L_tally, quant, L)
#         offsets = get_SA_offsets(SA_interval, red_SA, L, red_L_tally, index, alphabet, quant)
#         end_time = time.time()
#         runtimes_50_replicate.append(end_time-start_time)
        
#     runtimes_10.append(np.mean(runtimes_10_replicate))
#     runtimes_20.append(np.mean(runtimes_20_replicate))
#     runtimes_30.append(np.mean(runtimes_30_replicate))
#     runtimes_40.append(np.mean(runtimes_40_replicate))
#     runtimes_50.append(np.mean(runtimes_50_replicate))

# # plot running times:
# fig, ax = plt.subplots()
# sns.lineplot(x=ref_lengths, y=runtimes_50, ax=ax, label='read length = 50')
# sns.lineplot(x=ref_lengths, y=runtimes_40, ax=ax, label='read length = 40')
# sns.lineplot(x=ref_lengths, y=runtimes_30, ax=ax, label='read length = 30')
# sns.lineplot(x=ref_lengths, y=runtimes_20, ax=ax, label='read length = 20')
# sns.lineplot(x=ref_lengths, y=runtimes_10, ax=ax, label='read length = 10')
# plt.xlabel('ref length')
# plt.ylabel('runtime (s)')
# plt.tight_layout()
# plt.show()

################################################################
