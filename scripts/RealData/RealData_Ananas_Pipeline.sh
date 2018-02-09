#!/bin/bash

##USE: ./RealData_Pipeline.sh name extension input_fastq_directory output_fastq_directory Ananas_assembler_directory Ananas_results_directory Scripts_directory Overlap_value(OPTIONAL) ReadsDirection(OPTIONAL) n_parameter(OPTIONAL) n2_parameter(OPTIONAL) no_parameter(OPTIONAL) 

##EXAMPLE: ./RealData_Pipeline.sh my_name my_extension /updir/my_input_fastq_directory/ /updir/my_output_fastq_directory/ /updir/my_Ananas_assembler_directory/ /updir/my_Ananas_results_directory/ /updir/my_Scripts_directory/ 25(OPTIONAL) rf(OPTIONAL) 20(OPTIONAL) 2(OPTIONAL) 20(OPTIONAL)

#The file name of your fastq files should be name_1.fastq and name_2.fastq 

#The extension should be either fq or fastq

##The variable $name should specify the name of the fastq file you want to process. The name should be WITHOUT set of reads (_1/_2) and WITHOUT extension (.fastq or .fq) (for example: name1, name2).

##In the $in_fasta (directory with input files) there should be the fastq files. The fastq file names have the name of the sample, the set of reads (_1/_2) and the extension .fastq or .fq (e.g. name_1.fastq and name_2.fastq).


name=$1
extension=$2
in_fastq=$3
out_fastq=$4
Ananas=$5
Ananas_out=$6
scripts=$7

###optional#####
Ananas_ml=${8:-35}
Ananas_dir=${9:-fr}
Ananas_n=${10:-1}
Ananas_n2=${11:-1}
Ananas_no=${12:-1}


####prepare files####
sed "1~4 s/^.*/&\;${name}\/1/" ${in_fastq}${name}_1.${extension} > ${out_fastq}${name}_1_4Ananas.${extension}
sed "1~4 s/^.*/&\;${name}\/2/" ${in_fastq}${name}_2.${extension} > ${out_fastq}${name}_2_4Ananas.${extension}


################ANANAS###############

#####run Ananas#####
mkdir ${Ananas_out}

${Ananas}Ananas -i ${out_fastq}${name}_1_4Ananas.${extension},${out_fastq}${name}_2_4Ananas.${extension} -ml $Ananas_ml -dir $Ananas_dir -n $Ananas_n -n2 $Ananas_n2 -no $Ananas_no -outReadNames ReadName_Index_file_${name} -o ${Ananas_out} > ${Ananas_out}log.out

#####Get top from Ananas assembly####
${Ananas}GetTopFromFasta -if ${Ananas_out}final.fa -il ${Ananas_out}final.layout

