
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

for index, read in enumerate(SeqIO.parse("readsMappingToChr1.fa.txt", "fasta")):
    seq1 = SeqRecord(Seq(read.seq),
                     id="seq1")
    SeqIO.write(seq1, f"seq{index}.fasta", "fasta")
