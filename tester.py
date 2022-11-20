#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 14:32:47 2022

"""

import numpy as np
import subprocess
from io import StringIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIWWW
seq_iterator = SeqIO.parse("genomic.fna.txt", "fasta")
chromosome_1 = next(seq_iterator)
#chromosome 1 sequence
genome_seq = chromosome_1.seq
#print(genome_seq)
read_1 = next(SeqIO.parse("readsMappingToChr1.fa.txt", "fasta"))
read = SeqIO.parse("readsMappingToChr1.fa.txt", "fasta")
all_read = []
with open("example.fasta") as handle:
    for record in read:
        all_read.append(str(record.seq))
# Create two sequence files
seq1 = SeqRecord(Seq(read_1.seq),
                   id="seq1")
seq2 = SeqRecord(genome_seq,
                   id="seq2")
SeqIO.write(seq1, "seq1.fasta", "fasta")
SeqIO.write(seq2, "seq2.fasta", "fasta")

# Run BLAST and parse the output as XML
'''
result_handle = NCBIWWW.qblast("blastn", "nt", read_1.seq)
with open("my_blast.xml", "w") as out_handle:
 out_handle.write(result_handle.read())

result_handle.close()
result_handle = open("my_blast.xml")
blast_record = NCBIXML.read(result_handle)

cmd = "blastx -query opuntia.fasta -db nr -out opuntia.xml"
cmd += " -evalue 0.001 -outfmt 5"
subprocess.run(cmd, shell=True)
'''
output = NcbiblastnCommandline(query=seq1, subject=seq2,
                            evalue=0.001, remote=True, ungapped=True, outfmt=5)()[0]
blast_result_record = NCBIXML.read(StringIO(output))

for alignment in blast_result_record.alignments:
    for hsp in alignment.hsps:
        #if hsp.expect < 0.05:
            print("****Alignment****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value:", hsp.expect)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
