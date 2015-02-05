#!/usr/bin/perl -w
####Extract both not aligned paire-end reads
use strict;
use warnings;
use FileHandle;
use Getopt::Long;


my ($sam, $fasta, $output);
my $result = &GetOptions("sam|g=s{1}" => \$sam,
                         "fastq|f=s{1}" => \$fastq,
	            		"output|o=s{1}" => \$output);

unless ($result && defined($sam) && defined($fastq) && defined($output)) {
	print STDERR sprintf("perl %s --sam <sam alignments> --fasta <fasta file> --output <out>\n", $0);
	exit();
} # End of unless statement

# Reading alignments 
my %aligned;
my $fh = new FileHandle();
open ($fh, $sam) or die("Cannot open sam file\n");
while (<$fh>) {
	chomp;
	if($_=~/^([^@]\S+)\t/){ 	#sam file of d2b only contains aligned reads
		my $name=$1;
		$name=~s/\/[12]//; #remove the \1 or \2 of single-aligned reads
		$aligned{$name}=1;
	}
}
close ($fh);

# Reading fastq file
$fh = new FileHandle();
open ($fh, $fastq) or die("Cannot open fasta file\n");
my $name;

my $ofh = new FileHandle();
open ($ofh, sprintf(">%s", $output)) or die("Cannot create output fasta file\n");
my $count = 0;
my $record =0;
while (<$fh>) {
	chomp;
	if (length($_) != 0) {
		if ($_ =~ m/^@(MG\S+)(\/[12])$/) {
			$name=$1;
			if (defined($name) && !exists $aligned{$name}) {
			$count++;
			$record=1;
			print $ofh ">".$name."$2\n";
			} # End of if statement
		} # end of if statement
		elsif($record==1){
			print $ofh $_."\n";
			$record=0;
		}else{} # end of else statement
	} # end of if statement
} # End of while loop
close ($fh);
close ($ofh);
