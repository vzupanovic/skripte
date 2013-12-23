#!/bin/bash    
nucmer -maxmatch -c 100 -p nucmer genom.fasta contigi.fasta
show-coords -r -c -l nucmer.delta > nucmer.coords
show-snps -C nucmer.delta > nucmer.snps
show-tiling nucmer.delta > nucmer.tiling
delta-filter -m nucmer.delta > nucmer.delta.m
mummerplot nucmer.delta.m
