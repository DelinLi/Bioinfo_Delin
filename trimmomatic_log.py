#by Delin Li, delin.bio@gmail.com Jan 04, 2018
#Aim to summary trimmomatic log files
#python 3 ++ only

import argparse
import gzip
import sys
import re
import os

version = (3,0)
cur_version = sys.version_info
if(cur_version<version):
    sys.exit("require python 3.0 +, but yours is "+sys.version[:3]+"\n")



# Taking command line arguments from users
parser = argparse.ArgumentParser()
parser.add_argument('-l', '--logs',  help='log file directory', nargs='*', required=True)
parser.add_argument('-p', '--pair', help='suffix of log file',   default="p")
parser.add_argument('-o', '--output', help='output file', type=str, required=True)
args = parser.parse_args()


if os.path.exists(args.output) and os.path.getsize(args.output) > 0:
    fh = open(args.output, "a")
else:
    fh = open(args.output, "a")
    fh.write("File\tRaw_reads\tRaw_bp\tRaw_length\tTrimmed_reads\tTrimmed_reads_percen\tTrimmed_reads_pair\tTrimmed_bp\tTrimmed_bp_percen\tTrimmed_length\n")

if(args.pair=="p"):
    for file in args.logs:
        (raw_r_empty, raw_r, raw_l, raw_b, trimmed_r, trimmed_pr, trimmed_l, trimmed_b) = (0, 0, 0, 0, 0, 0, 0, 0)
        file_handler = 'false'
        if (file.endswith('.gz')):
            handle = gzip.open(file, 'rt')
            print("gz comrpressed file ")
            file_handler = "true"
        else:
            print("un-comrpressed file ")
            handle = open(file, 'r')
            file_handler = "true"
        while (file_handler):
            read1 = handle.readline().rstrip("\n\r")
            if len(read1) == 0:
                break
            else:
                array1 = list(map(int, read1.split(' ')[-4:]))
                read2 = handle.readline().rstrip("\n\r")
                array2 = list(map(int, read2.split(' ')[-4:]))
                if (array1[0] * array2[0] == 0):
                    if (array1[0] + array2[0] == 0):
                        raw_r_empty += 2
                    else:
                        raw_r_empty += 1
                        raw_r += 1
                        raw_b += (array1[2] + array1[3] + array2[2] + array2[3])
                        trimmed_b += (array1[0] + array2[0])
                        trimmed_r += 1
                else:
                    raw_r += 2
                    raw_b += (array1[2] + array1[3] + array2[2] + array2[3])
                    trimmed_b += (array1[0] + array2[0])
                    trimmed_r += 2
                    trimmed_pr += 2
        raw_l = raw_b / raw_r
        raw_r += raw_r_empty
        raw_b += raw_r_empty * raw_l
        if (file_handler):
            fh.write("%s\t%d\t%d\t%d\t%d\t%.3f\t%d\t%d\t%.3f\t%d\n" % (
            file, raw_r, raw_b, raw_l, trimmed_r, trimmed_r / raw_r, trimmed_pr, trimmed_b, trimmed_b / raw_b,
            trimmed_b / trimmed_r))
else: #signle end
    for file in args.logs:
        (raw_r_empty, raw_r, raw_l, raw_b, trimmed_r, trimmed_pr, trimmed_l, trimmed_b) = (0, 0, 0, 0, 0, 0, 0, 0)
        file_handler = 'false'
        if (file.endswith('.gz')):
            handle = gzip.open(file, 'rt')
            print("gz comrpressed file ")
            file_handler = "true"
        else:
            print("un-comrpressed file ")
            handle = open(file, 'r')
            file_handler = "true"
        while (file_handler):
            read1 = handle.readline().rstrip("\n\r")
            if len(read1) == 0:
                break
            else:
                array1 = list(map(int, read1.split(' ')[-4:]))
                if (array1[0]  == 0):
                        raw_r_empty += 1
                else:
                    raw_r += 1
                    raw_b += (array1[2] + array1[3])
                    trimmed_b += array1[0]
                    trimmed_r += 1
        raw_l = raw_b / raw_r
        raw_r += raw_r_empty
        raw_b += raw_r_empty * raw_l
        if (file_handler):
            fh.write("%s\t%d\t%d\t%d\t%d\t%.3f\t%d\t%d\t%.3f\t%d\n" % (
            file, raw_r, raw_b, raw_l, trimmed_r, trimmed_r / raw_r, trimmed_pr, trimmed_b, trimmed_b / raw_b,
            trimmed_b / trimmed_r))
fh.close