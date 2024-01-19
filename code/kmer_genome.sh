#!/bin/bash
echo "k,unique,distinct,total"
for ((k=1; k<=20; k++)); do
        #grep with -V flag so only process lines that doesn't start with "#"
        #awk command to find unique, distint and total number of kmers
        #unique kmer count from row that starts with 1 meaning X kmers repeat once
        #distinct kmer count from sum of 2nd column, total kmer from sum of distinct kmer times how many times it repeat
        grep -v '^#' ${k}mer_chr1.hist|awk '{
        if ($1==1) unique=$2;
        distinct+=$2;
        total+=$2*$1;
        }END{print "'$k',"unique","distinct","total}'
done
~             
