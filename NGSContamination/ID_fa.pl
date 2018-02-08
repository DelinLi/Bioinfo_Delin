#!/usr/bin/perl
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
while(<DATA>){
	if($_=~/^(\S+)\t(\S+)/){
		if($id_out{$1}){
			print"$1\t$2\n";
		}
	}
}
close DATA;



