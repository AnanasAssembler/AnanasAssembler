#!/bin/bash

##USE: ./Ananas_Pipeline.sh filename  input_fasta_directory output_fasta_directory scripts_directory nr_random_reads(for one pair member)

##EXAMPLE: ./Ananas_Pipeline.sh my_filename /updir/my_input_fasta_directory/ /updir/my_output_fasta_directory/ /updir/my_scripts_directory/ 1000(nr. of random reads I want for one member of the pair: here 1000 left reads, corresponding to 2000 [1000x2] reads in total)

#Dowload paired-end RNA-seq datasets (fasta format) from NCBI SRA repository using sratoolkit. Change file name to species1_1(2).fasta, species2_1(2).fasta, etc. Each file should specify a set of reads (left reads/_1 and right reads/_2).

##create a file ( my_filename ) in which you specify the names of the fasta files you want to process. The file names should be WITHOUT set of reads (_1/_2) and WITHOUT extension (.fasta) (for example: species1, species2), and each in a different line. Don't forget to leave an empty line at the end of the list file. This file should be in the same directory as the fasta files.

##In the $in_fasta (directory with input files) there should be the fasta files. The fasta file names have the name of the species, the set of reads (_1/_2) and the extension .fasta (e.g. species1_1.fasta and species1_2.fasta). The name of the fasta files (i.e. species1 and species1) will be use to correctly define the header of each sequence in the output fasta files (i.e. species1_1_4Ananas.fasta and species1_2_4Ananas.fasta)

###The Extract_ReadSubset_fromFasta.pl perl script should be in the script directory (/updir/my_script_directory/)


filename=$1
in_fasta=$2
out_fasta=$3
scripts=$4
nr_random_reads=$5

declare -a my_species_array
my_species_array=(`cat "$filename"`)
max=${#my_species_array[@]}

for ((i=0; i < $max; i++))
do
sed "s/>.*/&\;${my_species_array[$i]}\/1/" ${in_fasta}${my_species_array[$i]}_1.fasta | sed "/^>/s/ /\-/g" > ${out_fasta}${my_species_array[$i]}_1_4Ananas.fasta
sed "s/>.*/&\;${my_species_array[$i]}\/2/" ${in_fasta}${my_species_array[$i]}_2.fasta | sed "/^>/s/ /\-/g" > ${out_fasta}${my_species_array[$i]}_2_4Ananas.fasta
done

#####create the random subset with as many reads as you want coming from the same number of species you created the processable fasta files#######

####create the source files with all the reads###
cd $out_fasta
list_files_subset_1_4Ananas=$(find . -path '*1_4Ananas.fasta' | sort | paste -sd ' ')
list_files_subset_2_4Ananas=$(find . -path '*2_4Ananas.fasta' | sort | paste -sd ' ')

cat $list_files_subset_1_4Ananas > ALL_1.fasta
cat $list_files_subset_2_4Ananas > ALL_2.fasta

####get random reads and create corresponding files###

${scripts}Extract_ReadSubset_fromFasta.pl ${out_fasta}ALL_1.fasta ${out_fasta}ALL_2.fasta ${out_fasta}RandomSubset_1_4Ananas.fasta ${out_fasta}RandomSubset_2_4Ananas.fasta $nr_random_reads

sed -i "/^>/! s/ /\\n/g" ${out_fasta}RandomSubset_1_4Ananas.fasta
sed -i "/^>/! s/ /\\n/g" ${out_fasta}RandomSubset_2_4Ananas.fasta

rm ALL_1.fasta
rm ALL_2.fasta









