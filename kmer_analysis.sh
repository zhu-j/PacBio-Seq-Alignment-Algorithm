#!/bin/bash
mkdir -p /Users/jessicazhu/kmer_analysis
for ((i=1;i<=1054;i++)); do
	echo "k,unique,distinct,total" > /Users/jessicazhu/kmer_analysis/kmer_analysis_${i}.csv
	for ((k=1; k<=15; k++)); do
		#grep with -V flag so only process lines that doesn't start with "#"
		#awk command to find unique, distint and total number of kmers
		#unique kmer count from row that starts with 1 meaning X kmers repeat once
		#distinct kmer count from sum of 2nd column, total kmer from sum of distinct kmer times how many times it repeat
		grep -v '^#' /Users/jessicazhu/reads/${k}mer_${i}read.hist|awk '{ 
		if ($1==1) unique=$2; 
		distinct+=$2; 
		total+=$2*$1;
		}END{print "'$k',"unique","distinct","total}'
	done > /Users/jessicazhu/kmer_analysis/kmer_analysis_${i}.csv 
done
