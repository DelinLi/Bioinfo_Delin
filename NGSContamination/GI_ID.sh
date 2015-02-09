###BLAST Results
From GI to Taxid

for ACC in A00002 X53307 BB145968 CAA42669 V00181  AH002406  HQ844023
do
   echo -n -e "$ACC\t"
   curl -s "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACC}&rettype=fasta&retmode=xml" |\
   grep TSeq_taxid |\
   cut -d '>' -f 2 |\
   cut -d '<' -f 1 |\
   tr -d "\n"
   echo
 done


###Reason
$ curl "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=HQ844023&rettype=fasta&retmode=xml" 

http://www.ncbi.nlm.nih.gov/dtd/NCBI_TSeq.dtd">
<TSeqSet>
<TSeq>
  <TSeq_seqtype value="nucleotide"/>
  <TSeq_gi>341832806</TSeq_gi>
  <TSeq_accver>HQ844023.1</TSeq_accver>
  <TSeq_taxid>1056490</TSeq_taxid>
  <TSeq_orgname>Rotavirus A HC91xUK reassortant (UKg9KC-1)</TSeq_orgname>
  <TSeq_defline>Rotavirus A HC91xUK reassortant (UKg9KC-1) NSP3 protein gene, complete cds</TSeq_defline>
  <TSeq_length>942</TSeq_length>
  <TSeq_sequence>(...)</TSeq_sequence>
</TSeq>
</TSeqSet>
