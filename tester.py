#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 14:32:47 2022

@author: jessicazhu
"""

import numpy as np 
from Bio import SeqIO

seq_iterator = SeqIO.parse("genomic.fna.txt", "fasta")
chromosome_1 = next(seq_iterator)
#chromosome 1 sequence
genome_seq = chromosome_1.seq
print(genome_seq)
