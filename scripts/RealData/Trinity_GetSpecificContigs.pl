#!/usr/bin/perl

##use: perl_script.pl   Trinity_fasta     match_in_the_contig_header    Trinity_TOP.fasta(output)
###to get the first contigs, specify in the variable $match_contig as "c0_g1_i1"
###for example to get the "i1" contigs, specify "_i1"###

use warnings;
use strict; 

open(OUT, ">$ARGV[2]");
my $match_contig = $ARGV[1]; 

## 1. create a hash between the contig names and the sequences##
open(IN, $ARGV[0]);

my %contigs_sequence;
my $contig_name;
while (my $line = <IN>) {
    chomp($line);
    if ($line =~ m/^>(.*)$/) { $contig_name = $line; }
    else {
        chomp($line);
        push(@{$contigs_sequence{$contig_name}}, ($line));
    }
}

close IN;


## 2. Extract only the first contigs of each gene/isoform group###
for (grep /$match_contig/, keys %contigs_sequence) {
print OUT "$_\n@{$contigs_sequence{$_}}\n";
}

close OUT;

