#!/bin/bash
##USE: ./Statistics_Trinity.sh filename nr_species Reference_files_directory Trinity_assembly_directory Trinity_statistics_results_directory Scripts_directory BedTools2_directory

##EXAMPLE: ./Statistics_TrinityAndAnanas.sh /updir/my_filename TwoSpecies /updir/my_Reference_Genomes_directory/ /updir/my_Trinity_assembly_directory/ /updir/my_Trinity_statistics_results_directory/ /updir/my_scripts_directory/ /updir/BedTools2_directory/

#create a file ( my_filename ) in which you specify the names of the species included in the Trinity. Each name should be in a different line. Don't forget to leave an empty line at the end of the list file (same concept as the previous pipeline script). This file should be in the Reference files directory

#nr_species should specify how many species are included (for example TwoSpecies, ThreeSpecies, etc).

###The Reference_files_directory should contain:
#1) Reference Genome file for each species included in the analysis. Every reference genome should be names species_ReferenceGenome (same species name as in the my_filename). IMPORTANT: Replace the name of the assembly (generally a number in the header) with the name of the species.

######For each Reference Genome of the included species, you should download the corresponding gff file. Every gff file should be named species_gff.txt (same species name as in the my_filename). The gff file should be placed in the Reference files directory 

#####The reference genome and the corresponding gff file can be dowloaded from NCBI RefSeq

filename=$1
nr_species=$2
dir_reference=$3
dir_Trinity_assembly=$4
dir_Trinity_results=$5
dir_scripts=$6
dir_BedTools2=$7

####create the Trinity statitics results directories###
mkdir ${dir_Trinity_results}

######Get general statistics about the Trinity assembly#######
#####nr. of contigs####
grep '>' ${dir_Trinity_assembly}Trinity.fasta | wc -l > ${dir_Trinity_results}Trinity_number_of_contigs

#####bp assembled in the contigs######
${dir_scripts}Trinity_GetContigsLength.pl ${dir_Trinity_assembly}Trinity.fasta ${dir_Trinity_results}Trinity_contigs_length
awk "{ sum+=\$1} END {print sum}" ${dir_Trinity_results}Trinity_contigs_length > ${dir_Trinity_results}Trinity_TotalBpAssembled


####create the concatenated Reference Genome according to how many species are included in the analysis####
cd ${dir_reference}
list_ReferenceGenome=$(find . -path '*ReferenceGenome' | sort | paste -sd ' ')
cat $list_ReferenceGenome > ${dir_reference}${nr_species}_ReferenceGenome

####################
#######blast the Trinity assembly against the concatenated Reference Genome#####
makeblastdb -in ${dir_reference}${nr_species}_ReferenceGenome -parse_seqids -dbtype nucl

blastn -db ${dir_reference}${nr_species}_ReferenceGenome -query ${dir_Trinity_assembly}Trinity.fasta -num_threads 1 -outfmt 6 -out ${dir_Trinity_results}${nr_species}_Blast_Trinity_MinMatch


cd ${dir_reference}
declare -a my_species_array
my_species_array=(`cat "$filename"`)
max=${#my_species_array[@]}

for ((i=0; i < $max; i++))
do

##create the .bed files from the gff files###
grep -v "#" ${dir_reference}${my_species_array[$i]}_gff.txt | awk -v OFS="\t" "{print \$1,\$3,\$4}" | grep "^chr" > ${dir_reference}${my_species_array[$i]}_GeneAnnotation.bed


####Trinity####get Statistics based on Reference Genomes
grep "${my_species_array[$i]}" ${dir_Trinity_results}${nr_species}_Blast_Trinity_MinMatch > ${dir_Trinity_results}${nr_species}_Blast_Trinity_MinMatch_${my_species_array[$i]}

sort -k1,1 -k12,12nr -k11,11n ${dir_Trinity_results}${nr_species}_Blast_Trinity_MinMatch_${my_species_array[$i]} | sort -u -k1,1 --merge > ${dir_Trinity_results}${nr_species}_Blast_Trinity_MinMatch_${my_species_array[$i]}_best

sort -k9,9n -k10,10n ${dir_Trinity_results}${nr_species}_Blast_Trinity_MinMatch_${my_species_array[$i]}_best > ${dir_Trinity_results}${nr_species}_Blast_Trinity_MinMatch_${my_species_array[$i]}_best_sorted

${dir_scripts}Statistics_BlastDowstream.pl ${dir_Trinity_results}${nr_species}_Blast_Trinity_MinMatch_${my_species_array[$i]}_best_sorted ${dir_Trinity_results}${nr_species}_Blast_Trinity_MinMatch_${my_species_array[$i]}_best_sorted_OK

sort -k9,9n ${dir_Trinity_results}${nr_species}_Blast_Trinity_MinMatch_${my_species_array[$i]}_best_sorted_OK > ${dir_Trinity_results}${nr_species}_Blast_Trinity_MinMatch_${my_species_array[$i]}_start_file

${dir_scripts}Statistics_BlastGetBpGaps.pl ${dir_Trinity_results}${nr_species}_Blast_Trinity_MinMatch_${my_species_array[$i]}_start_file ${dir_Trinity_results}${my_species_array[$i]}_Trinity_BpGaps

awk -v OFS="\t" "{print \$9, \$10}" ${dir_Trinity_results}${nr_species}_Blast_Trinity_MinMatch_${my_species_array[$i]}_start_file > ${dir_Trinity_results}${nr_species}_Blast_Trinity_MinMatch_${my_species_array[$i]}_CoverageIntervals

perl -pe "s/^/chr\t/g" ${dir_Trinity_results}${nr_species}_Blast_Trinity_MinMatch_${my_species_array[$i]}_CoverageIntervals > ${dir_Trinity_results}${nr_species}_Blast_Trinity_MinMatch_${my_species_array[$i]}_CoverageIntervals.bed

sort -k1,1 -k2,2n ${dir_Trinity_results}${nr_species}_Blast_Trinity_MinMatch_${my_species_array[$i]}_CoverageIntervals.bed > ${dir_Trinity_results}${nr_species}_Blast_Trinity_MinMatch_${my_species_array[$i]}_CoverageIntervals_sorted.bed

${dir_BedTools2}mergeBed -i ${dir_Trinity_results}${nr_species}_Blast_Trinity_MinMatch_${my_species_array[$i]}_CoverageIntervals_sorted.bed > ${dir_Trinity_results}${nr_species}_Blast_Trinity_MinMatch_${my_species_array[$i]}_CoverageIntervals_sorted_NO_OVERLAPS.bed

${dir_BedTools2}intersectBed -a ${dir_Trinity_results}${nr_species}_Blast_Trinity_MinMatch_${my_species_array[$i]}_CoverageIntervals_sorted_NO_OVERLAPS.bed -b ${dir_reference}${my_species_array[$i]}_GeneAnnotation.bed -wo > ${dir_Trinity_results}${my_species_array[$i]}_Trinity_GeneIntersections

awk "{ sum+=\$7} END {print sum}" ${dir_Trinity_results}${my_species_array[$i]}_Trinity_GeneIntersections > ${dir_Trinity_results}${my_species_array[$i]}_Trinity_BpGenesOverlap

awk "{ sum+=\$4} END {print sum}" ${dir_Trinity_results}${nr_species}_Blast_Trinity_MinMatch_${my_species_array[$i]} > ${dir_Trinity_results}${my_species_array[$i]}_Trinity_TotalBpMatch

done

rm ${dir_Trinity_results}*_MinMatch*
rm ${dir_Trinity_results}*_GeneIntersections



