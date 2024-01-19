#!/bin/bash
for ((i=818; i<=1054; i++)); do
	echo "start ${i}"
	for ((k=1; k<=15; k++)); do
        	kat hist -o  /Users/jessicazhu/reads/${k}mer_${i}read.hist -m $k /Users/jessicazhu/comp561-project/seq${i}.fasta >/dev/null 2>&1
	done
	echo "end $i"
done
