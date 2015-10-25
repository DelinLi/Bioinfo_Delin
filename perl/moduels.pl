#!/use/bin/perl

###1. add comma to big number
####///////////////////////////////////////////////////////////
sub commify {
        my $text = reverse $_[0];
        $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
        return scalar reverse $text;
}
#//from online http://www.perlmonks.org/?node_id=110137
my $number = 987654321;
print STDOUT sprintf commify($number)."\n";
####///////////////////////////////////////////////////////////

###2. Check the out file exit or not
####///////////////////////////////////////////////////////////
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
####///////////////////////////////////////////////////////////


