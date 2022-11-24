# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 17:07:46 2022

@author: Connie
"""
from Bio import SeqIO
import numpy as np


'''
k_mer(read, k) returns a dictionary of k_mers of the 
input read with key being the index of k_mer and value being the k_mer
'''
def k_mer(read, k):
    size = len(read)-3+1
    output = {}
    for i in range(size-1):
        k_mer = read[i:k+i]
        output[i] = k_mer
    return output
        
'''
parser function to extract all reads from fasta file
and store into a list
'''
def parser(path):
    sequences = [] 
    for record in SeqIO.parse(path, "fasta"):
        read=str(record.seq)
        sequences.append(read)
    return sequences

def kmers_for_all_reads(reads):
    kmers = {}
    for index, read in enumerate(reads):
        kmers[index]=k_mer(read,3)
    return kmers

path='readsMappingToChr1.fa.txt'
R=parser(path)
S=kmers_for_all_reads(R)
