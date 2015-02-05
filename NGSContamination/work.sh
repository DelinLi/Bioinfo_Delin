###1. Extact the unaligned reads to fasta format
##The sam files only contains aligned reads
##Only un-aligned pair-end reads are exracted
perl extract_not_aligned.pl --sam <sam alignments> --fastq <fastq file> --output <out fasta file>

###2. Genrate a randomly dataset if the whole dataset is too large (like ~10M)
##if pair-end ,extract two respectively by specific number or percentage reads -> more randomly
perl extract_random_subset_fasta.pl --input <fastq files> --prefix <prefix string> --num <number/fraction>

###3. BLASTX and obtain the best hit
sh blast.sh <fasta file>
sort -k1,1 -k12,12nr -k11,11n Random.10k.P1_20000.fq.nt.blastn | sort -u -k1,1 --merge > Random.10k.P1_20000.Best

###4. From GI to Taxi_ID
awk 'BEGIN {FS="\t"} {print $2}' | sed 's/gi|\([0-9]*\)|.*/\1/' >Id.txt
perl ID_fa.pl -id Id.txt -data gi_taxid_nucl.dmp >IDs.taxi.txt
#gi_taxid_nucl.dmp from ftp://ftp.ncbi.nih.gov/pub/taxonomy/
awk 'BEGIN {FS="\t"} {print $2} IDs.taxi.txt' >IDs.taxi.only.txt
sort -u IDs.taxi.only.txt >IDs.taxi.only.unique.txt
##From taxi to organism http://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi

###5. Generate PieChart with R

