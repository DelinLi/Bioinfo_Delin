#!/bin/bash
#$ -M delin.bio@gmail.com
#$ -m ae
#$ -V 
#$ -l h_vmem=20G
#$ -pe smp 6
#$ -N ReadsCounts
#$ -cwd
#$ -l h_rt=80:00:00
#$ -o ReadsCounts.summary.out
#$ -e ReadsCounts.summary.err

featureCounts -T 6 -s 2 -a /datalus/guyongzhe/RNASeq_delin/genome/Gmax_275_Wm82.a2.v1.gene.gff3 -o Leaf.GmV2.CDS.readscount.20200328.txt  *bam   -g ID  -F GTF -t CDS -s 0 -p -B -C -f -O
