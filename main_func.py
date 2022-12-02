# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 16:13:44 2022
"""
from banded_needleman_wunsch import *
import _pickle as pickle
import pickle as pk
import ast
import time

def main(k, name):
    i = 0
    readPairs = {}
    for pair in kmerCountTable.keys():   
        read1 = R[pair[0]]
        read2 = R[pair[1]]
        L1 = len(read1)
        L2 = len(read2)
        common_kmer = kmerCountTable[pair]
        X = max((L1 - (common_kmer + k)), (L2 - (common_kmer + k)))
        gap_penalty = -1
        match_score = +1
        mismatch_score = -1
        print(i)
        R1R2_align = banded_needleman_wunsch(read1, read2, X, gap_penalty, match_score, mismatch_score)
        R2R1_align = banded_needleman_wunsch(read2, read1, X, gap_penalty, match_score, mismatch_score)
        # for overlap over 50% of the longer read sequence, we consider as predicted true read pair
        overlap_threshold = 0.5 * max(L1, L2)
        if max(R1R2_align, R2R1_align) > overlap_threshold:
            readPairs[pair] = (read1, read2)
        i+=1
        with open(name, 'wb') as txt_file:
            txt_file.write(pickle.dumps(readPairs))
        return readPairs
    
def evaluate(readPairs):    
    with open('true_pairs.pickle', 'rb') as handle:
        truePairs = pickle.load(handle)
    miss = 0
    correct = 0
    incorrect = 0
    total_readPairs = len(readPairs.keys())
    total_truePairs = len(truePairs.keys()) # truePairs = {(R1,R2): (R1 seq, R2 seq)}
    for pair in truePairs.keys():
        if pair in readPairs.keys() or (pair[1], pair[0]) in readPairs:
            correct +=1
    msg = ("Out of {} true pairs, {} were correct".format(total_truePairs, correct))
    print(msg)
    
    true = np.unique(list(sum(truePairs.keys())))
    for pair in readPairs.keys():
        if pair not in truePairs.keys() or (pair[1], pair[0]) not in truePairs.keys():
            if pair[0] in true or pair[1] in true:
                incorrect +=1
            else:
                miss +=1
                
    msg = ("Miss: {}, Incorrect: {}".format(miss, correct))
    return correct, incorrect, miss



# load list of sequencing reads
with open('reads_list.txt', 'rb') as handle:
    data = handle.read()
R = pk.loads(data)

# load common-kmer count table    
with open('kmerCountTable_threshold1500.txt', 'rb') as handle:
    data = handle.read()
kmerCountTable = pk.loads(data)
readPairs = {}

# time alignment 
start = time.time()
readPairs = main(k, name)
end = time.time()
alignment_time = end - start
print("Alignment duration: {}".format(alignment_time))
correct = incorrect = miss = 0
correct, incorrect, miss = evaluate(readPairs)