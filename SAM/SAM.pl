#!/usr/bin/perl

# Description: 
# Delin LI delin.bio@gmail.com

#V0.10 2014.10.17
use strict;
use warnings;
use FileHandle;
use Getopt::Long;

use constant VERSION => "0.10 (2015.04.1)";

use constant Homo => 4;
use constant Heter => 1;
use constant HeterPer => 0.2;
use constant HomoPer => 0.9;
use constant OverAllPer => 0.9;
use constant HeterTotal => 5;
use constant HeterCode =>3;
use constant RefCode =>1;
use constant AltCode =>2;
use constant MissingCode =>0;
use constant PreFix => "allele_counts";


my ($suffix, $countsDir, $outfile);
my ($homo,$heter,$homo_per,$heter_per,$overall_per,$heter_total,$heter_code,$ref_code,$alt_code,$missing_code);

my $result = &GetOptions("countsDir|f=s" =>\$countsDir,
			 "Postfix|pf=s"=> \$suffix,
			 "outfile|o=s"=> \$outfile,
			 "homo=i" => \$homo,
			 "heter=i" => \$heter,
			 "homo_per=i"=> \$homo_per,
			 "heter_per=i"=> \$heter_per,
			 "overall_per=i"=> \$overall_per,
			 "heter_total=i" => \$heter_total,
			 "heter_code=i" => \$heter_code,
			 "ref_code=s" => \$ref_code,
			 "alt_code=s" => \$alt_code,
			 "missing_code=s" => \$missing_code);

unless ($result && defined($suffix) && defined($countsDir) && defined($outfile)) {
	print STDERR sprintf("\n");
	print STDERR sprintf("perl %s --Postfix|-pf <Postfix of allele counts files> --countsDir|-f <counts dir> --outfile|-o <out file>\n", $0);
	print STDERR sprintf("\n");
	print STDERR sprintf("WHERE:\n");
	print STDERR sprintf("   --countsDir|-f <counts dir>    : Directory containing allele counts for each sample\n");
	print STDERR sprintf("   --Postfix|-pf  <counts dir>    : Postfix of allele counts files [DEFULT: %s]\n", PreFix);
	print STDERR sprintf("   --outfile|-o   <counts dir>    : output file name for genotyping\n");
	print STDERR sprintf("Optional\n");
	print STDERR sprintf("   --homo                         : min reads for a homozygous call [DEFAULT: %d]\n",Homo);
	print STDERR sprintf("   --heter                        : min reads of minor allele of a heterozygous call [DEFAULT: %d]\n",Heter);
	print STDERR sprintf("   --heter_total                  : min sum reads of two alleles of a heterozygous call [DEFAULT: %d]\n",HeterTotal);
	print STDERR sprintf("   --homo_per                     : min percentage for allele of a homozygous call [DEFAULT: %.1f]\n",HomoPer);
	print STDERR sprintf("   --heter_per                    : min percentage for alleles of a heterozygous call [DEFAULT: %.1f]\n",HeterPer);
	print STDERR sprintf("   --Overall_per                  : min percentage for all allele(s) of a heterozygous/homozygous call  [DEFAULT: %.1f]\n",OverAllPer);
	print STDERR sprintf("   --ref_code                     : Code of Reference allele [DEFAULT: %s]\n",RefCode);
	print STDERR sprintf("   --alt_code                     : Code of Alternative allele [DEFAULT: %s]\n",AltCode);
	print STDERR sprintf("   --missing_code                 : Code of Missiong [DEFAULT: %s]\n",MissingCode);
	print STDERR sprintf("\n");
	print STDERR sprintf("VERSION: %s\n", VERSION);
	print STDERR sprintf("\n");
	exit();
} # end of unless statement

$homo = Homo if (!defined($homo));        
$heter = Heter if (!defined($heter));                     
$homo_per = HomoPer if (!defined($homo_per));         
$heter_per = HeterPer if (!defined($heter_per));         
$overall_per = OverAllPer if (!defined($overall_per));         
$heter_total = HeterTotal if (!defined($heter_total));         
$ref_code = RefCode if (!defined($ref_code));                
$alt_code = AltCode if (!defined($alt_code));                
$heter_code = HeterCode if (!defined($heter_code));                
$missing_code = MissingCode if (!defined($missing_code));                
$suffix = PreFix if (!defined($suffix));                


##Check the out file exit or not
if(-e $outfile){
	my $answer;
	do{
		print STDERR sprintf("\n");
		print STDERR sprintf("WARNING: output '%s' already exists.\n",$outfile);
		print STDERR sprintf("\n");
		print STDERR sprintf("Would you like to delete the contents of that directory\n");
		print STDERR sprintf("and continue executing? (Y/N): ");
		$answer=<>;
		chomp($answer);
	} while ($answer !~ m/^(Y|N|yes|no)$/i);
	
	if ($answer =~ m/^(Y|yes)$/i) {
		print STDERR sprintf("\n");
		print STDERR sprintf(" o Removing existing '%s' ...\n",$outfile);
		my $command = sprintf("rm -f %s", $outfile);
		system($command);
		
		if ($? == -1) {
			print STDERR sprintf("FAILED\n");
			exit();
		}else{
			print STDERR sprintf("DONE\n");
		} 
	}else {
		print STDERR sprintf("\n");
		exit();
	}
}

print STDERR sprintf("\n");

