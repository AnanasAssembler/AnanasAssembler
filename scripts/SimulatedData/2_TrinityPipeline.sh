#!/bin/bash
##USE: ./Trinity_Pipeline.sh scripts_directory output_fasta_directory Trinity_assembler_directory Bowtie2_directory Trinity_output_directory Sequencing_type_fa_or_fq(OPTIONAL) Max_Memory_value(OPTIONAL) Trinity_CPUs_Value(OPTIONAL)

##EXAMPLE: ./Trinity_Pipeline.sh /updir/my_script_directory/ /updir/my_output_fasta_directory/ /updir/Trinity_assembler_directory/ /updir/Bowtie2_directory/ /updir/Trinity_output_directory/ fq(OPTIONAL) 20(OPTIONAL) 1(OPTIONAL)

##In the $out_fasta (directory with output fasta files) there should be the processable fasta files coming from the previous step.

###The Trinity_blast_downstream_lenient_4Pipeline.pl perl script should be in the script directory (/updir/my_script_directory/)

scripts=$1
out_fasta=$2
Trinity=$3
Bowtie=$4
Trinity_out=$5
Trinity_seqType=${6:-fa}
Trinity_maxMemory=${7:-40G}
Trinity_CPU=${8:-1}

cd $out_fasta
list_files_4Trinity_LEFT=$(find . -path '*_1_4Ananas.fasta' | sort | paste -sd,)
list_files_4Trinity_RIGHT=$(find . -path '*_2_4Ananas.fasta' | sort | paste -sd,)
echo $list_files_4Trinity_LEFT


PATH=$PATH:${bowtie}bowtie2
export PATH
export LC_CTYPE="en_US.UTF-8"

mkdir ${Trinity_out}

${Trinity}Trinity --seqType $Trinity_seqType --max_memory $Trinity_maxMemory --CPU $Trinity_CPU --no_normalize_reads --left $list_files_4Trinity_LEFT --right $list_files_4Trinity_RIGHT --output ${Trinity_out} > ${Trinity_out}log.out &&


#######run Trinity downstream######

cat ${out_fasta}*.fasta > ${Trinity_out}AllTogether.fasta

makeblastdb -in ${Trinity_out}Trinity.fasta -parse_seqids -dbtype nucl

blastn -db ${Trinity_out}Trinity.fasta -query ${Trinity_out}AllTogether.fasta -num_threads 5 -outfmt 6 -perc_identity 90 -out ${Trinity_out}Blast_back_Trinity &&

rm ${Trinity_out}AllTogether.fasta


${scripts}Trinity_blast_downstream.pl ${Trinity_out}Blast_back_Trinity ${Trinity_out}Trinity.fasta ${Trinity_out}Results_TrinityBlast_downstream_98 0.98 && 

sort -k1,1n -k2,2nr -k4,4nr ${Trinity_out}Results_TrinityBlast_downstream_98 | awk -F"\t" '!_[$1]++' > ${Trinity_out}Results_TrinityBlast_downstream_98_SORTED &&


######R scripts######
R --vanilla <<EOF

input <- read.table("${Trinity_out}Results_TrinityBlast_downstream_98_SORTED", sep="\t", head=F)

############histogram contig length####
input\$log_contig_length <- log(input\$V2)
min_value_log <- min(input\$log_contig_length)
max_value_log <- max(input\$log_contig_length)
pdf("${Trinity_out}histogram_contig_length_log.pdf")
hist(input\$log_contig_length, xlab = "log contig length", col = "red", breaks = 50, xlim = c(min_value_log, max_value_log), main="Trinity")
dev.off()

###calculate max ratio###
input\$max_ratio <- input\$V4/input\$V5

max_ratio_eq1 <- length(which(input\$max_ratio == 1))
max_ratio_all <- length(input\$max_ratio)
accuracy <- max_ratio_eq1/max_ratio_all

write.table(max_ratio_eq1, "${Trinity_out}max_ratio_eq1.txt")
write.table(max_ratio_all, "${Trinity_out}max_ratio_all.txt")
write.table(accuracy, "${Trinity_out}accuracy.txt") 


#Correlation between contig length and max ratio#
species_col <- as.factor(input\$V3)
pdf("${Trinity_out}corr_contig_length_max_ratio_diff_spec.pdf")
plot(input\$max_ratio~input\$V2, xlab = "contig length", ylab = "max ratio", col=species_col, pch=16, main="Trinity")
legend("bottomright", legend=levels(species_col), pch=16, col=seq_along(levels(species_col)))
dev.off()

#Correlation between COVERAGE (nr. reads) and max ratio#
species_col <- as.factor(input\$V3)
pdf("${Trinity_out}corr_coverage_max_ratio_diff_spec.pdf")
plot(input\$max_ratio~input\$V5, xlab = "coverage", ylab = "max ratio", col=species_col, pch=16, main = "Trinity")
legend("bottomright", legend=levels(species_col), pch=16, col=seq_along(levels(species_col)))
dev.off()


#Correlation between Contig length and coverage (nr. reads)#
species_col <- as.factor(input\$V3)
pdf("${Trinity_out}corr_contig_length_coverage_diff_spec.pdf")
plot(input\$V5~input\$V2, xlab = "contig length", ylab = "coverage", col=species_col, pch=16, main = "Trinity")
legend("bottomright", legend=levels(species_col), pch=16, col=seq_along(levels(species_col)), cex = 0.9)
dev.off()

EOF


