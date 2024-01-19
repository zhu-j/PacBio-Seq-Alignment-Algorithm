#!/bin/bash
cd /Users/jessicazhu/kmer_analysis
for file in $(echo *.csv); do
    gnuplot <<EOF
    set terminal pngcairo color
    set output "/Users/jessicazhu/kmer_png_plot/output_${file}.png"
    set key bottom right
    set xlabel "k"
    set ylabel "kmer count"
    set datafile separator ','  #csv file
    plot "$file" u 1:2 t 'unique' w lines linewidth 3 linecolor "red", "$file" u 1:3 t 'distinct' w lines linewidth 3 linecolor "green", "$file" u 1:4 t 'total' w lines linewidth 3 linecolor "blue" 
EOF
done
