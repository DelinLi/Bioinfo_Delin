#!/usr/bin/perl 
##By Delin LI  delin.bio@gamil.com
###Extract Gene/mRNA/CDS sequence from reference genome to fasta format
##V1.00 only work on specific chromosome/scaffold for now
##V1.01 could extract sequence from any chromosomes, the input gene table should be with 4 columns

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
	print STDERR sprintf("\n");
	print STDERR sprintf("perl %s --reference|-f <reference fasta> --genes|-g <gene table with four columns> --output|-o <out fasta>\n", $0);
	print STDERR sprintf("\n");
	print STDERR sprintf("\t--genes|-g The four columns of gene table: Chr, ID, StartPosition, EndPosition\n\n");
	exit();
} 

my (%genome,$seq,$name);
my $fh = new FileHandle();
open ($fh, $reference) or die("Cannot open fasta file\n");
while (<$fh>) {
	chomp;
	if (length($_) != 0) {
		if ($_ =~ m/^>(\S+)/){
            if(defined($name)){
                $genome{$name}=$seq;
            }
        $name=$1;
        $seq="";
        }else{
            $seq.=$_;
        }
    }#end of if
}#end of while
close ($fh);

my $ofh = new FileHandle();
open (IN,$genes) or die("Cannot open gene table\n");
open ($ofh,">$output") or die("Cannot create output fasta file\n");
while(<IN>){
	my@line=split(/\t/,$_);
	if(!defined($genome{$line[0]})){
		print STDERR "Line without recognized chromosome name, please ignore this if this is a headline:\n";
		print STDERR "$_";
	}else{
		my$sequence=substr($genome{$line[0]}, $line[2] - 1, $line[3] - $line[2] + 1);
		print $ofh ">$line[1]\n";
		print $ofh "$sequence\n";
	}
}#end of while loop
close IN;
close $ofh;
