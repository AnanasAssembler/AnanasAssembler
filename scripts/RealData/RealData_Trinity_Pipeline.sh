#!/bin/bash

##USE: ./RealData_Pipeline.sh name extension input_fastq_directory output_fastq_directory  Trinity_assembler_directory Bowtie2_directory Trinity_results_directory Scripts_directory Sequencing_type_fa_or_fq(OPTIONAL) Max_Memory_value(OPTIONAL) Trinity_CPUs_Value(OPTIONAL)

##EXAMPLE: ./RealData_Pipeline.sh my_name my_extension /updir/my_input_fastq_directory/ /updir/my_output_fastq_directory/ /updir/my_Trinity_assembler_directory/ /updir/my_Bowtie2_directory/ /updir/my_Trinity_Results_directory/ /updir/my_Scripts_directory/ fa(OPTIONAL) 20G(OPTIONAL) 1(OPTIONAL)

#The file name of your fastq files should be name_1.fastq and name_2.fastq 

#The extension should be either fq or fastq

##The variable $name should specify the name of the fastq file you want to process. The name should be WITHOUT set of reads (_1/_2) and WITHOUT extension (.fastq or .fq) (for example: name1, name2).

##In the $in_fasta (directory with input files) there should be the fastq files. The fastq file names have the name of the sample, the set of reads (_1/_2) and the extension .fastq or .fq (e.g. name_1.fastq and name_2.fastq).


name=$1
extension=$2
in_fastq=$3
out_fastq=$4
Trinity=$5
Bowtie=$6
Trinity_out=$7
scripts=$8

###optional#####
Trinity_seqType=${9:-fq}
Trinity_maxMemory=${10:-40G}
Trinity_CPU=${11:-1}


####prepare files####
sed "1~4 s/^.*/&\;${name}\/1/" ${in_fastq}${name}_1.${extension} > ${out_fastq}${name}_1_4Ananas.${extension}
sed "1~4 s/^.*/&\;${name}\/2/" ${in_fastq}${name}_2.${extension} > ${out_fastq}${name}_2_4Ananas.${extension}

#####run Trinity######

PATH=$PATH:${bowtie}bowtie2
export PATH
export LC_CTYPE="en_US.UTF-8"

mkdir ${Trinity_out}

${Trinity}Trinity --seqType $Trinity_seqType --max_memory $Trinity_maxMemory --CPU $Trinity_CPU --no_normalize_reads --left ${out_fastq}${name}_1_4Ananas.${extension} --right ${out_fastq}${name}_2_4Ananas.${extension} --output ${Trinity_out} > ${Trinity_out}log.out &&

####Get top from Trinity assembly###

${scripts}Trinity_GetSpecificContigs.pl ${Trinity_out}Trinity.fasta c0_g1_i1 ${Trinity_out}Trinity_TOP.fasta

sed -i "/^>/! s/ /\\n/g" ${Trinity_out}Trinity_TOP.fasta







