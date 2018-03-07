#by Delin Li, delin.bio@gmail.com Feb 04, 2018
#Aim to extract fasta sequence
#python 3 ++ only
#updated Feb 08, 2018 
from Bio import SeqIO
import argparse
import sys
import re

version = (3,0)
cur_version = sys.version_info
if(cur_version<version):
    sys.exit("require python 3.0 +, but yours is "+sys.version[:3]+"\n")

# Taking command line arguments from users
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fasta', help='fasta file path', type=str, required=True)
parser.add_argument('-l', '--list', help='four column list input file to extract', type=str, required=False)
parser.add_argument('-o', '--output', help='output file', type=str, required=True )
args = parser.parse_args()


fw=open(args.output,'w') #open the output file
#read in list
listF=open(args.list,'rU')
ids=dict() #creat a dict to store the list: id: information(start end )
chrs={}
for row in listF:
    rows=row.split("\t")
    ids[rows[0]]={}
    ids[rows[0]]["chr"] = rows[1]
    ids[rows[0]]["start"] = rows[2]
    ids[rows[0]]["end"] = rows[3]
    chrs[rows[1]]=1

listF.close()
# read in the fasta file and extract & output them
inFasta=open(args.fasta,'rU') #open the fasta file
for record in SeqIO.parse(inFasta,'fasta'):
    id_tem=re.sub(r'\s.*$','',record.id)
    if id_tem in chrs:
        print(id_tem)
        for key,_ in ids.items():
            if(ids[key]["chr"] == id_tem):
                fw.write(">" + key + "\n")
                start=int(ids[key]["start"])-1
                end=int(ids[key]["end"])
                fw.write(str(record.seq)[start:end]+"\n")

fw.close() # close output file handle
