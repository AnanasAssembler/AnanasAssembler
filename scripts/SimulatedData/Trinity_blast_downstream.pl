#!/usr/bin/perl

##use: perl_script.pl input(BLAST output) input2(Trinity fasta file) out_file percentage_dec 

##output format: name of the Trinity contig \t length of the contig \t name of the species in the contig \t how many reads from that species are in the contig \t number of total reads in the contig##

use warnings;
use strict; 

my $percentage_dec = $ARGV[3];
open(OUT, ">$ARGV[2]");



##1. Create a hash between the Trinity contig names and how many species map there##

open(IN, $ARGV[0]);

my %TriBla_contig_species;

while (my $line = <IN>) {
    chomp($line);
    my @species = ($line =~ m/;(.*?)\/[1-2]\t/);
    my @length = ($line =~ m/length=(.*?);/);
    my @col = split("\t", $line); 
        for my $elem (@length) {
            if (($col[3] le $elem) && ($col[3] >= ($elem * $percentage_dec)) and ($col[2] <= (($col[3]/$col[3])*100.00) && ($col[2] >= (($col[3]/$col[3])*100.00) * $percentage_dec))) {
                push(@{ $TriBla_contig_species{$col[1]} }, @species);
        }
    }
}

close(IN);


##2. Create a hash between the Trinity contig names and their lengthe##
open(IN_2, $ARGV[1]);

my %TriBla_contig_length;
while (my $line = <IN_2>) {
    chomp($line);
    if ($line =~ m/^>/) {    
        my @split = split(" ", $line);
          foreach (@split) {
            s/>//;
            s/>//;
          }
       my @length = ($line =~ m/len\=(.*?)\s/);
       $TriBla_contig_length{$split[0]} = "@length";  
    }
}

close(IN_2);


##3. Count the number of the different species for each contig in the hash ###
for my $contig (keys %TriBla_contig_species) {
    my $total_count;
    my %count;     
        for my $species (@{$TriBla_contig_species{$contig}}) {
        $total_count++;
        $count{$species}++;
    }
    foreach (keys %count) {
         print OUT "$contig\t$TriBla_contig_length{$contig}\t$_\t$count{$_}\t$total_count\n";
    }
}

close(OUT);

