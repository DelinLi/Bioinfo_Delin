#!/usr/bin/perl -w
##By Delin LI  delin.bio@gamil.com
###Extract Gene/mRNA/CDS sequence from reference genome to fasta format
##V1.00 only work on specific chromosome/scaffold for now

use strict;
use warnings;
use FileHandle;
use Getopt::Long;

use constant CONTEXT_LENGTH => 100;	
use constant true => 1;
use constant false => 0;

my ($output,  $reference,$genes);
my $result = &GetOptions("reference|f=s{1}" => \$reference,
                         "genes|g=s{1}" => \$genes,
	            		"output|o=s{1}" => \$output);

#Ask for parameters
unless ($result && defined($reference) && defined($genes) && defined($output)) {
	print STDERR sprintf("perl %s --reference|f <reference fasta> --genes|g <gene table with three columns> --output <out fasta>\n", $0);
	exit();
} 


my $seq;
my $fh = new FileHandle();
open ($fh, $reference) or die("Cannot open fasta file\n");
my $target=0;
while (<$fh>) {
	chomp;
	if (length($_) != 0) {
		if ($_ =~ m/^>DAUCA_BGI_Solexa_pseudochr_2_1_Chr_08/){
			$target=1;
		}elsif($_ =~ m/^>/){ $target=0}
		else{
			if($target==1){
				$seq .= $_;
			} # End of else statement
		}
	} # end of if statement
} # end of while loop
close ($fh);

my $ofh = new FileHandle();
open (IN,$genes) or die("Cannot open gene table\n");
open ($ofh,">$output") or die("Cannot create output fasta file\n");
while(<IN>){
	my@line=split(/\t/,$_);
	my$sequence=substr($seq, $line[1] - 1, $line[2] - $line[1] + 1);
	print $ofh ">$line[0]\n";
	print $ofh "$sequence\n";
}
close IN;
close $ofh;
