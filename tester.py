#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 14:32:47 2022

"""
"""
Prequisite: install blast+ suite from https://anaconda.org/bioconda/blast
"""
import numpy as np
import subprocess
from io import StringIO
from Bio.Blast import NCBIXML
from itertools import islice, count
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIWWW
import _pickle as pickle
seq_iterator = SeqIO.parse("genomic.fna.txt", "fasta")
chromosome_1 = next(seq_iterator)

#chromosome 1 sequence
genome_seq = chromosome_1.seq
chr1_seq = SeqRecord(genome_seq,
                 id="chr1_seq")
SeqIO.write(chr1_seq, "chr1_seq.fasta", "fasta")
dict_reads = {}

#blast individual reads against chromosome 1 subject seq
#generate dictionary where key has format (index, read sequence) and value of (genome region start, genome region end)
for index,read in enumerate(SeqIO.parse("readsMappingToChr1.fa.txt", "fasta")):
    seq1 = SeqRecord(Seq(read.seq),
                     id="seq1")
    SeqIO.write(seq1, "seq1.fasta", "fasta")
    (seq_s_e,blank) = NcbiblastnCommandline(query="seq1.fasta", subject="seq2.fasta", outfmt='10 sstart send')()
    list_s_e = seq_s_e.split("\n")
    pair_s_e = list_s_e[0]
    start_end = pair_s_e.split(",")
    dict_reads[(index,str(read.seq))] = start_end
    break
with open('test_case.txt', 'wb') as txt_file:
    txt_file.write(pickle.dumps(dict_reads))

#find pairs of read that corresponds to same genomic region
#code doesn't work right now
dict_altered = dict_reads
read_i = ""
read_pairs = []
for key,value in dict_reads.items():
    search_region = dict_altered.pop(key)
    read_i = key[1]
    continue
    if (value == search_region):
        (i,j) = (read_i, key[1])
        read_pairs.append(i,j)
