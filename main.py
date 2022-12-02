# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 16:13:44 2022
"""
from k_mer import *
from banded_needleman_wunsch import *
import _pickle as pickle
import pickle as pk
import ast

# 1: Generate a table containing kmers (k=10) for all reads
path = "readsMappingToChr1.fa.txt"
R = parser(path)
k = 10
kmerTable = kmers_for_all_reads(R, k)
kmerDict, kmerDictWithIndex = kmerDict(kmerTable)

# 2: Generate a common-kmers count table for all read pairs
s = 1500
kmerCountTable = kmerFreqPerPair(kmerDict, s)

# 3 For all read pairs with common-kmer counts > some threshold, we perform alignment
readPairs = {}
for pair in kmerCountTable.keys():
    read1 = R[pair[0]]
    read2 = R[pair[1]]
    L1 = len(read1)
    L2 = len(read2)
    common_kmer = kmerCountTable[pair]
    X = (L1 - (common_kmer + k)) + (L2 - (common_kmer + k))  # maximum possible of gaps
    gap_penalty = -1
    match_score = +1
    mismatch_score = -1
    R1R2_align = banded_needleman_wunsch(read1, read2, X, gap_penalty, match_score, mismatch_score)
    R2R1_align = banded_needleman_wunsch(read2, read1, X, gap_penalty, match_score, mismatch_score)
    # for overlap over 50% of the longer read sequence, we consider as this pair as coming from same region in genome
    overlap_threshold = 0.5 * max(L1, L2)
    if max(R1R2_align, R2R1_align) > overlap_threshold:
        readPairs[pair] = (read1, read2)

readID_pairs = list(readPairs.keys())
readSeq_pairs = list(readPairs.values())
print("sequencing reads of the same origin: ", readSeq_pairs)

# store the read pairs found by alignment
with open('read_pairs.txt', 'wb') as file:
    file.write(pickle.dumps(readPairs))

# evaluate how many read pairs were correctly found
with open("true_pairs.txt", "r") as data:
    truePairs = ast.literal_eval(data.read())

miss = 0
correct = 0
incorrect = 0
total_readPairs = len(readPairs.keys())
total_truePairs = len(truePairs.keys())  # truePairs = {(R1,R2): (R1 seq, R2 seq)}

for pair in truePairs.keys():
    if pair in readPairs.keys() or (pair[1], pair[0]) in readPairs:
        correct += 1
msg = ("Out of {} true pairs, {} were correct".format(total_truePairs, correct))
print(msg)

true = np.unique(list(sum(truePairs.keys())))
for pair in readPairs.keys():
    if pair not in truePairs.keys() or (pair[1], pair[0]) not in truePairs.keys():
        if pair[0] in true or pair[1] in true:
            incorrect += 1
        else:
            miss += 1





















