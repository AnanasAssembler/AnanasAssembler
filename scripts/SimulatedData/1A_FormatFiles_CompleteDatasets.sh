#!/bin/bash

##USE: ./Ananas_Pipeline.sh filename input_fasta_directory output_fasta_directory

##EXAMPLE: ./Ananas_Pipeline.sh my_filename /updir/my_input_fasta_directory/ /updir/my_output_fasta_directory/

#Dowload paired-end RNA-seq datasets (fasta format) from NCBI SRA repository using sratoolkit. Change file name to species1_1(2).fasta, species2_1(2).fasta, etc. Each file should specify a set of reads (left reads/_1 and right reads/_2).

##create a file ( my_filename ) in which you specify the names of the fasta files you want to process. The file names should be WITHOUT set of reads (_1/_2) and WITHOUT extension (.fasta) (for example: species1, species2), and each in a different line. Don't forget to leave an empty line at the end of the list file. This file should be in the same directory as the fasta files.

##In the $in_fasta (directory with input files) there should be the fasta files. The fasta file names have the name of the species, the set of reads (_1/_2) and the extension .fasta (e.g. species1_1.fasta and species1_2.fasta). The name of the fasta files (i.e. species1 and species1) will be use to correctly define the header of each sequence in the output fasta files (i.e. species1_1_4Ananas.fasta and species1_2_4Ananas.fasta)


filename=$1
in_fasta=$2
out_fasta=$3

declare -a my_species_array
my_species_array=(`cat "$filename"`)
max=${#my_species_array[@]}

for ((i=0; i < $max; i++))
do
sed "s/>.*/&\;${my_species_array[$i]}\/1/" ${in_fasta}${my_species_array[$i]}_1.fasta | sed "/^>/s/ /\-/g" > ${out_fasta}${my_species_array[$i]}_1_4Ananas.fasta
sed "s/>.*/&\;${my_species_array[$i]}\/2/" ${in_fasta}${my_species_array[$i]}_2.fasta | sed "/^>/s/ /\-/g" > ${out_fasta}${my_species_array[$i]}_2_4Ananas.fasta
done

