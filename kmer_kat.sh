#!/bin/bash
echo "k,unique,distinct,total"
for ((k=1; k<=20; k++)); do
        kat hist -o ${k}mer_chr1.hist -m $k /Users/jessicazhu/comp561-project/chr1_seq.fasta
done
~    
