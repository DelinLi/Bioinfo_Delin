#!/bin/bash
#$ -M delin.bio@gmail.com
#$ -m ae
#$ -V 
#$ -l h_vmem=24G
#$ -pe smp 10
#$ -N Trimming
#$ -wd /datalus/guyongzhe/RNASeq_delin/alignment
#$ -t 1-652:1
#$ -tc 20
#$ -l h_rt=60:00:00
#$ -o $JOB_NAME_$TASK_ID.out
#$ -e $JOB_NAME_$TASK_ID.err

proc=10
Ref="/datalus/guyongzhe/RNASeq_delin/Gmax_275_v2.db"
db="GmV2"
Prefixs=($(find *leaf.trimmed-paired-1.fq.gz | awk '{ sub(/.trimmed-paired-1\.fq\.gz$/, "", $0); print $0; }' | sort -u))

prefix="${Prefixs[$SGE_TASK_ID - 1]}"

if [ ! -e "$prefix"."$db".Aligned.sortedByCoord.out.bam ]
then
	echo $prefix
	STAR --runThreadN 12 --genomeDir ../Gmax_275_v2.db/  --readFilesIn "$prefix".trimmed-paired-1.fq.gz "$prefix".trimmed-paired-2.fq.gz \
	--outFilterMultimapNmax 10  --alignIntronMax 50000  --outFileNamePrefix "$prefix"."$db". --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat
fi
