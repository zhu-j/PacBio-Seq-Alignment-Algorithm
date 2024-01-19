from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
import _pickle as pickle

seq_iterator = SeqIO.parse("genomic.fna.txt", "fasta")
chromosome_1 = next(seq_iterator)

genome_seq = chromosome_1.seq
chr1_seq = SeqRecord(genome_seq,
                     id="chr1_seq")
SeqIO.write(chr1_seq, "chr1_seq.fasta", "fasta")
dict_reads = {}

# blast individual reads against chromosome 1 subject seq
# generate dictionary where key has format (index, read sequence) and value of (mismatch, gaps)
for index, read in enumerate(SeqIO.parse("readsMappingToChr1.fa.txt", "fasta")):
    seq1 = SeqRecord(Seq(read.seq),
                     id="seq1")
    SeqIO.write(seq1, "seq1.fasta", "fasta")
    (mismatch_gap, blank) = NcbiblastnCommandline(query="seq1.fasta", subject="seq2.fasta", outfmt='10 mismatch gaps')()
    list_mismatch_gap = mismatch_gap.split("\n")
    pair_mismatch_gap = list_mismatch_gap[0]
    errors = pair_mismatch_gap.split(",")
    dict_reads[(index, str(read.seq))] = errors

with open('errors.txt', 'wb') as txt_file:
    txt_file.write(pickle.dumps(dict_reads))
