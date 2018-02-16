#!/usr/bin/perl
##use: perl_script.pl TWO_column_layout_top_file_GLOB length_of_the_contigs

use warnings;
use strict; 
use Data::Dumper;

open(OUT, ">$ARGV[1]");

## 1. create a hash between the contig names and the consensus read indexes (new index, NO raw)##
open(IN, $ARGV[0]);

my %contigs_length;
my $contig_name;

while (my $line = <IN>) {
    if ($line =~ m/^>(.*)$/) {
        $contig_name = $1;
    }
    else {
            chomp($line);
  	    my @all = split(" ", $line);
            $contigs_length{$contig_name} = $all[1];  #get the last value of the arrays;it does it automatically because of the loop
    }
}

#print the contigs length values##
for my $contig (keys %contigs_length) {
    print OUT "$contigs_length{$contig}\n";
}

close IN;
close OUT;

