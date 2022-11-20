# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 17:07:46 2022

@author: Connie
"""
from Bio import SeqIO
import numpy as np


def k_mer(sequence, k):
    size = len(sequence)-3+1
    output = []
    for i in range(size-1):
        k_mer = sequence[i:k+i]
        output.append(k_mer)
        print(i, k_mer)
    return output
        
def parser():
    sequences = [] 
    for record in SeqIO.parse("readsMappingToChr1.fa.txt", "fasta"):
        sequences.append(record)


s="GTAGAGCTGT "
out=k_mer(s, 3)