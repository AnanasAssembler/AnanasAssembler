#!/usr/bin/perl

##use: perl_script.pl TWO_column_layout_top_file_GLOB ReadName_Index_file outfile_name
##output format: name of the contig \t length of the contig \t name of the species in the contig \t how many reads from that species are in the contig \t number of total reads in the contig##

use warnings;
use strict; 


open(OUT, ">$ARGV[2]");


## 1. create a hash between the contig names and the consensus read indexes (new index, NO raw)##
open(IN, $ARGV[0]);

my %contigs_readIndex;
my %contigs_length;
my $contig_name;
while (my $line = <IN>) {
    if ($line =~ m/^>(.*)$/) { $contig_name = $1; }
    else {
        chomp($line);
	my @all = split(" ", $line);
        push(@{$contigs_readIndex{$contig_name}}, ($all[0]));
        $contigs_length{$contig_name} = $all[1];    #get the last value of the arrays;it does it automatically because of the loop
    }
}


close IN;


## 2. create a hash between the consensus read indexes and the number of species (how many raw reads collapsed in the consensus reads) ##
open(IN_READ_names, $ARGV[1]);

my $string;
my %readIndex_species;
while (my $line = <IN_READ_names>) {                                                         
    chomp $line;   
    my @all = split(" ", $line);
    foreach (@all) {
      s/\(//;
      s/\)//;
    }
        if ($line =~ m/ \- (.*)$/) {$string = $1;}    
        my @species = ($string =~ m/[0-9]\;(.*?)\//g);
        $readIndex_species{$all[1]} = [@species];
}				

close IN_READ_names;


## 3. create a hash between the contig names and the number of species (basically the species replace the consensus read index ###
my %contigs_species;
for my $contigName (keys %contigs_readIndex) {
    for my $consensusId (@{$contigs_readIndex{$contigName}}) {     
        if(exists $readIndex_species{$consensusId}) {
            push(@{$contigs_species{$contigName}}, @{$readIndex_species{$consensusId}});
        }
    }	
}


## 4. count the number of the different species for each contig in the hash ###
for my $contig_last (keys %contigs_species) {      
    my $total_count; #count of all reads(species) collapsed in the consensus reads
    my %count; #hash relating each species name and its count 
    for my $species_last (@{$contigs_species{$contig_last}}) {
        $total_count++;
        $count{$species_last}++;
    }
    foreach (keys %count) {
        print OUT "$contig_last\t$contigs_length{$contig_last}\t$_\t$count{$_}\t$total_count\n";
    }
}

close OUT;

