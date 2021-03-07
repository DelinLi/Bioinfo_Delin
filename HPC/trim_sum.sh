#!/bin/bash
#$ -M delin.bio@gmail.com
#$ -m ae
#$ -V 
#$ -l h_vmem=2G
#$ -pe smp 1
#$ -N Trimming
#$ -cwd
#$ -l h_rt=80:00:00
#$ -o Trimming.summary.out
#$ -e Trimming.summary.err

python trim.log.py  -l *log -o trimming.summary.20200326.txt
