from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIWWW
import _pickle as pickle


# blast individual reads against chromosome 1 subject seq
# generate dictionary where key has format (index, read sequence) and value of (genome region start, genome region end)
read = next(SeqIO.parse("readsMappingToChr1.fa.txt", "fasta"))
seq1 = SeqRecord(Seq(read.seq),
                 id="read1")
SeqIO.write(seq1, "read1.fasta", "fasta")
