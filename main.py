# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 16:13:44 2022
"""
from k_mer import *
from banded_needleman_wunsch import *
import _pickle as pickle
import pickle as pk
import ast
import time

#1: load list of sequencing reads
with open('reads_list.txt', 'rb') as handle:
    data = handle.read()
R = pk.loads(data)

#2: load common-kmer table 
with open('kmerCountTable_threshold1500.txt', 'rb') as handle:
    data = handle.read()
kmerCountTable = pk.loads(data)

start = time.time()
#3: For all read pairs with common-kmer counts > some threshold, we perform alignment
readPairs = {}
counter = 0
for pair in kmerCountTable.keys():
    read1 = R[pair[0]]
    read2 = R[pair[1]]
    L1 = len(read1)
    L2 = len(read2)
    common_kmer = kmerCountTable[pair]
    k = 10
    X = (L1 - (common_kmer + k)) + (L2 - (common_kmer + k))  # maximum possible of gaps
    gap_penalty = -1
    match_score = +1
    mismatch_score = -1
    counter+=1
    R1R2_align = banded_needleman_wunsch(read1, read2, X, gap_penalty, match_score, mismatch_score)
    R2R1_align = banded_needleman_wunsch(read2, read1, X, gap_penalty, match_score, mismatch_score)
    # for overlap over 50% of the longer read sequence, we consider as this pair as coming from same region in genome
    overlap_threshold = 0.5 * min(L1, L2)
    if max(R1R2_align, R2R1_align) > overlap_threshold:
        readPairs[pair] = (read1, read2)
        #print(pair)
    #print("one pass")

readID_pairs = list(readPairs.keys())
readSeq_pairs = list(readPairs.values())
#print("sequencing reads of the same origin: ", readSeq_pairs)
end = time.time()
duration = end - start

start_file = time.time()
with open('read_p.txt', 'wb') as output:
    pk.dump(readPairs, output, protocol=pk.HIGHEST_PROTOCOL)
end_file = time.time()
file_duration = end_file - start_file
with open('duration2.txt', 'w') as f:
    f.write('running_time = ' + str(duration) + '\n')
    f.write('start = ' + str(start) + '\n')
    f.write('end = ' + str(end) + '\n')
    f.write('file_write_time = ' + str(file_duration) + '\n')
    f.write('start of file write = ' + str(start_file) + '\n')
    f.write('end of file write= ' + str(end_file) + '\n')

#print(duration)
# store the read pairs found by alignment
with open('read_id.txt', 'wb') as out_file:
    pk.dump(readID_pairs, out_file, protocol=pk.HIGHEST_PROTOCOL)

with open('read_seq_value.txt', 'wb') as o:
    pk.dump(readID_pairs, o, protocol=pk.HIGHEST_PROTOCOL)

with open('read_pairs.txt', 'wb') as file:
    file.write(pk.dumps(readPairs, protocol=pk.HIGHEST_PROTOCOL))

with open('read_same_origin.txt', 'wb') as file:
    file.write(pk.dumps(readSeq_pairs, protocol=pk.HIGHEST_PROTOCOL))


"""
# evaluate how many read pairs were correctly found
with open("true_pairs.txt", "rb") as data:
    truePairs = ast.literal_eval(pk.load(data))
    truePairs = pk.load(data)
    #truePairs =
    #print(truePairs)

miss = 0
correct = 0
incorrect = 0
total_readPairs = len(readPairs.keys())
total_truePairs = len(truePairs.keys())  # truePairs = {(R1,R2): (R1 seq, R2 seq)}

true = np.unique(list(sum(truePairs.keys())))
for pair in readPairs.keys():
    if pair in truePairs.key() or (pair[1], pair[0]) in truePairs.keys():
        correct += 1
    else:
       if pair not in truePairs.keys() or (pair[1], pair[0]) not in truePairs.keys():
            if pair[0] in true or pair[1] in true:
                 incorrect += 1
            else:
                 miss += 1

with open('sensitivity.txt', 'w') as sens:
    sens.write('miss = ' + str(miss) + '\n')
    sens.write('incorrect = ' + str(incorrect) + '\n')
    sens.write('correct = ' + str(correct) + '\n')
    sens.write('total pairs = ' + str(total_readPairs) + '\n')
    sens.write('total true pairs = ' + str(total_truePairs) + '\n')
    sens.write(str(msg) + '\n')"""





















