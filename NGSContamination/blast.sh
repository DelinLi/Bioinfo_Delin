#!/bin/bash

dbpath="PathofNT" 
db="nt"
###NCBI’s NR/NT (Non-Redundant) database, Dowload from ftp://ftp.ncbi.nlm.nih.gov/blast/db/
proc=1
hits=5
evalue="1e-10"

if [ $# -eq 0 ]
then
	echo ""
	echo "ERROR: Please specify the path of sequence file for alignment"
	echo ""
	exit 0;
else
	input=$1
	output="$input"."$db".blastn

	time -v -o "$input"."$db".blastn.resources \
        blastn -db "$dbpath"/"$db" -query "$input" -out "$output" -evalue "$evalue" -show_gis -num_threads "$proc" -outfmt "6 qseqid sseqid pident qlen length mismatch gapope evalue bitscore" 
        # -max_target_seqs 1     --- https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/bty833/5106166?redirectedFrom=fulltext
	#			-num_alignments "$hits" -num_descriptions "$hits" -num_threads "$proc" -outfmt 6 # -max_target_seqs 1
fi
