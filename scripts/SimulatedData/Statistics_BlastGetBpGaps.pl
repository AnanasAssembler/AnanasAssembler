#!/usr/bin/perl

##use: perl_script.pl SpeciesName_start_file output_file

use warnings;
	use strict; 

open(IN, "$ARGV[0]");
open(OUT, ">$ARGV[1]");

my @coordinates;
my @coordinates1;
my $max=0; ###start from 0####
my $gap=0;
my $last_value; #####last coordinate;

while ( my $line = <IN>) {
    chomp($line);
    my @all = split("\t", $line);
    push(@coordinates, $all[8]);
    push(@coordinates1, $all[9]);
    $last_value = $all[9];
}         


foreach my $coord ( 0..$#coordinates ) {
    my $start = $coordinates[$coord];
    my $end = $coordinates1[$coord];
        if ($max == 0) {
            $gap = $start;
            $max=$end;
        } else {
            if($start <= $max && $end >= $max) {
                $max=$end;
            }
            if ($start > $max && $end > $max) { 
                    $gap+=($start - $max);
                    $max=$end;
            }
        }
}
print OUT "Maximum_value_reached\tBp_gaps_reached\n$max\t$gap\n";

