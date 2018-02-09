#!/usr/bin/perl

##use: perl_script.pl Trinity_or_Ananas_vs_REF_best_sorted output_file

use warnings;
use strict; 

open(OUT, ">>$ARGV[1]");


## 1. create a hash between the contig names and the consensus read indexes (new index, NO raw)##
open(IN, $ARGV[0]);


while (my $line = <IN>) {
    chomp $line;
    my @all = split("\t", $line);
    if ($all[8] > $all[9]) {
    print OUT "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]\t$all[6]\t$all[7]\t$all[9]\t$all[8]\t$all[10]\t$all[11]\n";
    }    
    if ($all[8] <= $all[9]) {
    print OUT "$all[0]\t$all[1]\t$all[2]\t$all[3]\t$all[4]\t$all[5]\t$all[6]\t$all[7]\t$all[8]\t$all[9]\t$all[10]\t$all[11]\n";
    }
}
close IN;
close OUT;
