#!/usr/bin/perl

#Used for notify the process via email
#Nov 5, 2014
#Delin Li 

#####There is a bug in the module line  /System/Library/Perl/Extras/5.16/IO/Socket/SSL.pm It's a bug!!!
#"m{^(!?)(?:(SSL(?:v2|v3|v23|v2/3))|(TLSv1[12]?))}i" from "m{^(!?)(?:(SSL(?:v2|v3|v23|v2/3))|(TLSv1[12]?))$}i"

use strict;
use warnings;
use Net::SMTP::TLS; ##install libauthen-sasl-perl
use Getopt::Long;

use constant FROM => "delin.email\@gmail.com";
use constant OUTGOING_SERVER => "smtp.gmail.com";
use constant USER => "delin.email\@gmail.com";
use constant PASS => "delin\@data";

#use constant FROM => "delin_email\@sina.com";
#use constant OUTGOING_SERVER => "smtp.sina.com";
#use constant USER => "delin_email\@sina.com";
#use constant PASS => "delin\@data2bio";

#use constant FROM => "delin_email\@163.com";
#use constant OUTGOING_SERVER => "smtp.163.com";
#use constant USER => "delin_email\@163.com";
#use constant PASS => "delin\@data2bio";

use constant RECIPIENT => "delin.bio\@data2bio.com";

my $template;
my($log,$email,$job);
my $result = &GetOptions("log|f=s" =>\$log,
                         "email|e=s"=> \$email,
                         "job|j=s" => \$job);

$email = RECIPIENT if (!defined($email));    

if(defined($log)){
	open IN,"<$log" or die "$!\n";
	while(<IN>){
        	$template.=$_;
	}
	close IN;
}else{
	$template="";
} 

my $server = Net::SMTP::TLS->new( OUTGOING_SERVER,
                                                                  Hello => OUTGOING_SERVER,
                                                                  Port => 587, #25,465 or 587
                                                                  User => USER,
                                                                  Password => PASS
                                                                 );  

$server->mail(FROM);
$server->to($email);
$server->data();
$server->datasend(sprintf("To: %s\n", $email));
$server->datasend(sprintf("From: %s\n", FROM));
$server->datasend(sprintf("Subject: Job %s is done\n", $job));
$server->datasend("\n");
$server->datasend(sprintf("%s\n", $template));
$server->dataend;
$server->quit;

sleep(2);       # Sleep for 2 seconds


