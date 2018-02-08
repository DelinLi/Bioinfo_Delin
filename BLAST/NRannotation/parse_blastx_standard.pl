#!/usr/bin/perl -w

use strict;
use warnings;
use FileHandle;
use Getopt::Long;

use constant true => 1;
use constant false => 0;
use constant TOTAL_HITS => 5;

my ($inFile, $outFile, $hits);
my $result = &GetOptions("input|i=s{1}" => \$inFile, "output|o=s{1}" => \$outFile, "hits|h:i{1}" => \$hits);

unless ($result && defined($inFile) && defined($outFile)) {
	print STDERR sprintf("\n");
	print STDERR sprintf("perl %s --input <blast file> --output <output>\n", $0);
	print STDERR sprintf("\n");
	exit();
} # end of unless statement

$hits = TOTAL_HITS if (!defined($hits) || $hits !~ m/^\d+$/);

my $fh = new FileHandle();
open ($fh, $inFile) or die("Cannot open input file\n");

my $ofh = new FileHandle();
open ($ofh, sprintf(">%s", $outFile)) or die("Cannot create output file\n");

print $ofh sprintf("# %s\n", scalar(localtime(time)));
print $ofh sprintf("# Input: %s\n", $inFile);
print $ofh sprintf("# Hits Count: %s\n", $hits);
print $ofh sprintf("#\n");
print $ofh sprintf("# No.\tQuery\tQuery Length");

for (my $i=0; $i < $hits; $i++) {
	print $ofh sprintf("\tHit Description #%s\tE-Value #%s", $i + 1, $i + 1);
} # end of for loop
print $ofh sprintf("\n");

my $index = 0;
my $query;
my $queryLength;
my $desc;
my $eval;
my $count;
my %hits;
while (<$fh>) {
	chomp;
	if (length($_) != 0) {
		if ($_ =~ m/^Query= /) {	# New Entry
			if (defined($count) && $count != 0) {	# Need to print
				print $ofh sprintf("%s\t%s\t%s", ++$index, $query, $queryLength);
				for (my $i=0; $i < $hits; $i++) {
					if (exists $hits{$i + 1}) {
						print $ofh sprintf("\t%s\t%s", $hits{$i + 1}->{"desc"}, $hits{$i + 1}->{"eval"});
					} # end of if statement
					else {
						print $ofh sprintf("\t--\t--");
					} # end of else statement
				} # end of for loop
				print $ofh sprintf("\n");
			} # end of if statemnet

			$query = $';
			undef($queryLength);
			undef($desc);
			undef($eval);
			$count = 0;
			%hits = ();	# Reset hits hash

			# Continue to read until "Length=" is found
			my $stop = false;
			while (!$stop) {
				my $line = <$fh>;	chomp($line);
				if ($line =~ m/^Length=(\d+)$/i) {
					$queryLength = $1;
					$stop = true;
				} # end of if statement
			} # end of while loop
		} # end of if statement
		elsif ($_ =~ m/^\*{1,} No hits found \*{1,}$/i) {	# No hits found
			print $ofh sprintf("%s\t%s\t%s\t%s\n", ++$index, $query, $queryLength, $&);
		} # end of if statement
		elsif ($_ =~ m/^>/) {	# Found a hit description
			$desc = $';
			
			# increment count
			$count++;

			# Continue to read until "Length=" is found
			my $stop = false;
			while (!$stop) {
				my $line = <$fh>;	chomp($line);
				if ($line =~ m/^Length=(\d+)/) {
					$stop = true;
				} # end of if statement
				elsif (length($line) != 0) {	# Append description
					$desc .= sprintf(" --%s", $line);
				} # end of if statement
			} # end of while loop

			# Save description
			$hits{$count}->{"desc"} = $desc;
		} # end of else if statement
		elsif ($_ =~ m/^\s+Score =/ && $_ =~ m/Expect(\S+)? = (\S+)/) {
			$eval = $2;
			
			# Save eval
			$hits{$count}->{"eval"} = $eval;
		} # end of else if statement
	}  # End of if statement
} # end of while loop

# Last record
if (defined($count) && $count != 0) { # need to print
	print $ofh sprintf("%s\t%s\t%s", ++$index, $query, $queryLength);
	for (my $i=0; $i < $hits; $i++) {
		if (exists $hits{$i + 1}) {
			print $ofh sprintf("\t%s\t%s", $hits{$i + 1}->{"desc"}, $hits{$i + 1}->{"eval"});
		} # end of if statement
		else {
			print $ofh sprintf("\t--\t--");
		} # end of else statement
	} # end of for loop
	print $ofh sprintf("\n");
} # end of if statemnet

close ($fh);
close ($ofh);
