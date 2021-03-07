#!/bin/bash
#$ -M delin.bio@gmail.com
#$ -m ae
#$ -V 
#$ -l h_vmem=4G
#$ -pe smp 10
#$ -N Trimming
#$ -wd /datalus/guyongzhe/RNASeq_delin/trimmed
#$ -t 1-652:1
#$ -tc 20
#$ -l h_rt=2:00:00
#$ -o $JOB_NAME_$TASK_ID.out
#$ -e $JOB_NAME_$TASK_ID.err

WORK_DIR="/datalus/guyongzhe/RNASeq_delin/trimmed"
proc=10
jar="/datalus/guyongzhe/opt/Trimmomatic-0.39_DL/trimmomatic-0.39.jar"
Prefixs=($(find *_1_clean.fq.gz | awk '{ sub(/_1_clean\.fq\.gz$/, "", $0); print $0; }' | sort -u))

prefix="${Prefixs[$SGE_TASK_ID - 1]}"

if [ ! -e "$prefix".leaf.trimmed-paired-1.fq.gz ]
then
	echo $prefix
         java -jar $jar PE -threads $proc -phred33 -trimlog ${prefix}.leaf.trimmed.log "$prefix"_1_clean.fq.gz "$prefix"_2_clean.fq.gz "$prefix".leaf.trimmed-paired-1.fq.gz\
         "$prefix".leaf.trimmed-singletons-1.fq.gz  "$prefix".leaf.trimmed-paired-2.fq.gz "$prefix".leaf.trimmed-singletons-2.fq.gz \
         LEADING:15 TRAILING:15 SLIDINGWINDOW:4:20 MINLEN:100 TOPHRED33
fi
