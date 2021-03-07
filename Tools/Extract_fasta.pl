#!/usr/bin/perl
###Extract fasta sequece from targeted "data" based on "id" files with one ID per line
### Delin LI
###Apr 2, 2015

use warnings;
use strict;
use Getopt::Long;

my $id_file="";
my $results;
my $data;
$results=GetOptions("id=s" =>\$id_file,"data=s"=> \$data);
my %id_out;
open ID,"<$id_file" or die "$!";
while(<ID>){
	if($_=~/^(\S+)/){
		$id_out{$1}=1;
	}else{
		print"warnings";
	}
}
close ID;

open DATA,"<$data" or die"$!";
my$status;
while(<DATA>){
	if($_=~/^>(\S+)/){ ##be careful about the space in this line
		if($id_out{$1}){
			$status=1;
			print"$_";
		}else{
			$status=0;
		}
	}else{
		if($status==1){
			print"$_";
		}
	}
}
close DATA;



