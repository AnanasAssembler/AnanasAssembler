#!/bin/bash

##USE: ./Ananas_Pipeline.sh scripts_directory output_fasta_directory Ananas_assembler_directory Ananas_output_directory Overlap_value(OPTIONAL) ReadsDirection(OPTIONAL) n_parameter(OPTIONAL) n2_parameter(OPTIONAL) no_parameter(OPTIONAL)

##EXAMPLE: ./Ananas_Pipeline.sh /updir/my_script_directory/ /updir/my_output_fasta_directory/ /updir/Ananas_assembler_directory/ /updir/Ananas_output_directory/ 25(OPTIONAL) rf(OPTIONAL) 15(OPTIONAL) 2(OPTIONAL) 15(OPTIONAL)

##In the $out_fasta (directory with output fasta files) there should be the processable fasta files coming from the previous step.

###The Ananas_dowstream.pl perl script should be in the script directory (/updir/my_script_directory/)

scripts=$1
out_fasta=$2
Ananas=$3
Ananas_out=$4
Ananas_ml=${5:-35}
Ananas_dir=${6:-fr}
Ananas_n=${7:-1}
Ananas_n2=${8:-1}
Ananas_no=${9:-1}

mkdir ${Ananas_out}

cd $out_fasta
list_files_4Ananas=$(find . -path '*.fasta' | sort | paste -sd,) 

${Ananas}Ananas -i $list_files_4Ananas -ml $Ananas_ml -dir $Ananas_dir -n $Ananas_n -n2 $Ananas_n2 -no $Ananas_no -outReadNames ReadName_Index_file -o ${Ananas_out} > ${Ananas_out}log.out &&


${Ananas}GetTopFromFasta -if ${Ananas_out}final.fa -il ${Ananas_out}final.layout &&


#######run Ananas downstream######


grep -v -e '<SCAFFOLD>' -e '<CONTIG_READCOUNT>' -e '<CONTIG_PAIRCOUNT>' -e '<\/CONTIG>' -e '<SCAFFOLD_READCOUNT>' -e '<SCAFFOLD_PAIRCOUNT>' -e '\/SCAFFOLD\>' -e '<unknown>' ${Ananas_out}final.layout.top | sed '/^<CONTIG>/{s/\t/-/g}' | awk '{print $1, $5}' | sed '/^</{s/^<CONTIG>-//g}' | sed '/^Contig/{s/-.*-.*$//g}' > ${Ananas_out}TWO_column_layout_top_file_GLOB &&


${scripts}Ananas_dowstream.pl ${Ananas_out}TWO_column_layout_top_file_GLOB ${Ananas_out}ReadName_Index_file ${Ananas_out}Results_Ananas_downstream_top && 

sort -k1,1n -k2,2nr -k4,4nr ${Ananas_out}Results_Ananas_downstream_top | awk -F"\t" '!_[$1]++' > ${Ananas_out}Results_Ananas_downstream_top_SORTED &&


######R scripts######
R --vanilla <<EOF

input <- read.table("${Ananas_out}Results_Ananas_downstream_top_SORTED", sep="\t", head=F)

############histogram contig length####
input\$log_contig_length <- log(input\$V2)
min_value_log <- min(input\$log_contig_length)
max_value_log <- max(input\$log_contig_length)
pdf("${Ananas_out}histogram_contig_length_log.pdf")
hist(input\$log_contig_length, xlab = "log contig length", col = "red", breaks = 50, xlim = c(min_value_log, max_value_log), main="Ananas top 35 bp overlap")
dev.off()

###calculate max ratio###
input\$max_ratio <- input\$V4/input\$V5

max_ratio_eq1 <- length(which(input\$max_ratio == 1))
max_ratio_all <- length(input\$max_ratio)
accuracy <- max_ratio_eq1/max_ratio_all

write.table(max_ratio_eq1, "${Ananas_out}max_ratio_eq1.txt")
write.table(max_ratio_all, "${Ananas_out}max_ratio_all.txt")
write.table(accuracy, "${Ananas_out}accuracy.txt") 


#Correlation between contig length and max ratio#
species_col <- as.factor(input\$V3)
pdf("${Ananas_out}corr_contig_length_max_ratio_diff_spec.pdf")
plot(input\$max_ratio~input\$V2, xlab = "contig length", ylab = "max ratio", col=species_col, pch=16, main="Ananas top 35 bp overlap")
legend("bottomright", legend=levels(species_col), pch=16, col=seq_along(levels(species_col)))
dev.off()

#Correlation between COVERAGE (nr. reads) and max ratio#
species_col <- as.factor(input\$V3)
pdf("${Ananas_out}corr_coverage_max_ratio_diff_spec.pdf")
plot(input\$max_ratio~input\$V5, xlab = "coverage", ylab = "max ratio", col=species_col, pch=16, main = "Ananas top 35 bp overlap")
legend("bottomright", legend=levels(species_col), pch=16, col=seq_along(levels(species_col)))
dev.off()


#Correlation between Contig length and coverage (nr. reads)#
species_col <- as.factor(input\$V3)
pdf("${Ananas_out}corr_contig_length_coverage_diff_spec.pdf")
plot(input\$V5~input\$V2, xlab = "contig length", ylab = "coverage", col=species_col, pch=16, main = "Ananas top 35 bp overlap")
legend("bottomright", legend=levels(species_col), pch=16, col=seq_along(levels(species_col)), cex = 0.9)
dev.off()

EOF


