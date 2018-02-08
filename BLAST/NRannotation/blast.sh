#!/bin/bash

dbpath="NTDatabase"
db="nt"
proc=20
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
                 blastn -db "$dbpath"/"$db" -query "$input" -out "$output" -evalue "$evalue" -show_gis \
                                -num_alignments "$hits" -num_descriptions "$hits" -num_threads "$proc"
fi


dbpath="/home/shared/references/blast/"
db="nr"
proc=20
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
        output="$input"."$db".blastx

        time -v -o "$input"."$db".blastx.resources \
                 blastx -db "$dbpath"/"$db" -query "$input" -out "$output" -evalue "$evalue" -show_gis \
                                -num_alignments "$hits" -num_descriptions "$hits" -num_threads "$proc"
fi