#!/usr/bin/perl
use strict;
use warnings;

use Bio::DB::Fasta;
use FileHandle;
use Getopt::Long;

use constant VERSION => "1.0 (2015.12.17)";

#use constant FILTER => 1;

my ($fasta, $outfile,$region);

my $result = &GetOptions("fasta|f=s" =>\$fasta,
			"region|r=s" =>\$region,
			 "out|o=s" => \$outfile);

unless (($result && defined($fasta) && defined($outfile) && defined($region)) || scalar(@ARGV)==4 ) {
	print STDERR sprintf("\n");
	print STDERR sprintf("perl %s --fasta|-f <fasta file>  --region|-r <regions> --outfile|-o <out file>\n", $0);
	print STDERR sprintf("\nor\n");
	print STDERR sprintf("perl %s <fasta file> <Chr> <Start> <End>\n",$0);
	print STDERR sprintf("WHERE:\n");
	exit();
}

#$fasta = FILTER if (!defined($fasta));

if(!($result && defined($fasta) && defined($outfile) && defined($region)) and scalar(@ARGV)==4){
	my($ref,$chr,$start,$end)=@ARGV;
	my $db = Bio::DB::Fasta->new($ref);
	my $seq = $db->seq($chr, $start, $end);
	print $seq,"\n";
	exit();
}


#my $db = Bio::DB::Fasta->new($file, -makeid => \&makeid);

#sub makeid {
#	my ($head) = @_;
#	$head =~ /^>([^:]+)/ or die qq(Invalid header "$head");
#		$1;
#}

my%Geno;
open IN,"<$region" or die "$!\n";
open OUT,">>$outfile" or die "$!\n";
my$count_snps=0;
my $db = Bio::DB::Fasta->new($fasta);
while(<IN>){
	chomp($_);
	my($chr,$start,$end)=split(/\t/,$_);
	my $seq = $db->seq($chr, $start, $end);
	if(defined($seq)){
		print OUT $chr,"\t",$start,"\t",$seq,"\n";
		$count_snps ++;
	}
}
close IN;
close OUT;
print "Total sites with reference allele: $count_snps\n";
