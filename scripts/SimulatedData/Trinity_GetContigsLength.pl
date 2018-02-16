#!/usr/bin/perl

##use: perl_script.pl input(Trinity fasta file) out_file  

use warnings;
use strict; 
#use Data::Dumper;

open(OUT, ">$ARGV[1]");

open(IN, $ARGV[0]);


my %Tri_contigs_length;
while (my $line = <IN>) {
    chomp($line);
    if ($line =~ m/^>/) {    
        my @split = split(" ", $line);
          foreach (@split) {
            s/>//;
            s/>//;
          }
       my @length = ($line =~ m/len\=(.*?)\s/);
       $Tri_contigs_length{$split[0]} = "@length";  
    }
}

#print the contigs length values##
for my $contig (keys %Tri_contigs_length) {
    print OUT "$Tri_contigs_length{$contig}\n";
}

close(IN);
close(OUT);

