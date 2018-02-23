#!/bin/bash
#$ -S /bin/bash
#$ -N kallisto
#$ -o $HOME/run/logs
#$ -e $HOME/run/logs
module load hisat2/2.0.5
module load samtools/1.3.1
module load stringtie/1.3.1c
##The 1.3.3b version of stringtie is installed in home directory and located in $HOME/bin
module load picard  
module load kallisto/0.43.1
module load igenome-human/hg38
## Usage
kallisto_index="/mnt/work1/users/bhklab/Users/zhaleh/hg38/kallisto_hg38.idx"

## CCLE: $1 study, $2 name of cell line, $3 bam file
## GNE: $1 study, $2 name of cell line, $3 input_fastq_1, $second_fastq input_fastq_2
## UHN: $1 study, $2 name of sample, $3 input_fastq_1, $second_fastq input_fastq_2
## GRAY: $1 study, $2 name of sample, $3 /mnt/work1/users/bhklab/Data/Gray/rnaseq/SRP026537/$sample_name/$sample_name.sra
## GEO/GSE_id: $1 study, $2 SRR_id
if [ $# -eq 5 ]; then
  study_name=$1
  sample_name=$2
  first_fastq=$3
  second_fastq=$4
  output_dir=$5
  alignment_flag=$6 #for kallisto alignment flag set to FALSE else set to TRUE
  alignment_tool=$7
  quant_tool=$8
else
  study_name=$1
  sample_name=$2
  input_file=$3
  output_dir=$4
  alignment_flag=$5 #for kallisto alignment flag set to FALSE else set to TRUE
  alignment_tool=$6
  quant_tool=$7
fi


# Navigate to file
if [ ! -d $output_dir ]; then
  mkdir  $output_dir
fi

cd  $5
if [ ! -d "$study_name/$sample_name" ]; then
	mkdir -p $study_name/$sample_name
fi
date

if [ ${input_file: -4} == ".bam" ]; then
	echo "Running picard"
  echo $study_name
  echo $sample_name
  echo $inpute_file
  java -Xmx16g -jar $picard_dir/picard.jar SamToFastq  I=$inpute_file FASTQ=$study_name/$sample_name/$sample_name_1.fastq SECOND_END_FASTQ=$study_name/$sample_name/$sample_name_2.fastq
  first_fastq=$study_name/$sample_name/$sample_name_1.fastq
  second_fastq=$study_name/$sample_name/$sample_name_2.fastq
  fastq_rm=TRUE
fi

if [ ${input_file: -4} == ".sra" ]; then
  echo"Running fastq-dump"
  echo $study_name
  echo $sample_name
  echo $inpute_file

	fastq-dump --split-files $inpute_file --outdir $study_name/$sample_name/
  first_fastq=$study_name/$sample_name/$sample_name_1.fastq
  second_fastq=$study_name/$sample_name/$sample_name_2.fastq
  fastq_rm=TRUE
fi

date
if $alignment_flag; then
  echo "Running kallisto"
  echo $study_name
  echo $sample_name
  echo $first_fastq
  echo $second_fastq
  kallisto quant -t 8 -i  -o $study_name/$sample_name $first_fastq $second_fastq
else if [$alignment_tool == "HISAT"]; then
  echo "Running HISAT2"
  hisat2 -p 12 --dta -x /mnt/work1/users/bhklab/Users/zhaleh/Genome/GRCh38/Hisat/grch38_tran/grch38_tran -1 $first_fastq -2 $second_fastq -S $study_name/$sample_name/Aligned.out.sam
  date
  echo "Sort and Convert Sam file to bam"
  samtools sort -@ 8 -o $study_name/$sample_name/Aligned.out.sorted.bam $study_name/$sample_name/Aligned.out.sam
  echo "Running stringtie"
  date
  stringtie $study_name/$sample_name/Aligned.out.sorted.bam -v -o $study_name/$sample_name/stringtie_output.gtf -A $study_name/$sample_name/gene_abund.tab -p 8 -G /mnt/work1/users/bhklab/Users/zhaleh/Genome/GRCh38/Gencode/gencode.v26.annotation.gtf
fi

if $fastq_rm; then
	rm $first_fastq
	rm $second_fastq
fi
date

#stringtie $study_name/$sample_name/Aligned.out.sorted.bam -v -o $study_name/$sample_name/ballgown/stringtie_output.gtf -e -B -p 8 -G $HOME/Genome/GRCh38/Gencode/gencode.v26.annotation.gtf
#stringtie $study_name/$sample_name/Aligned.out.sorted.bam -v -o $study_name/$sample_name/test/stringtie_output.gtf -p 8 -G $HOME/Genome/GRCh38/Gencode/gencode.v26.annotation.gtf
#date

#echo "Running Rail-rna"
#module load igenome-human/hg38
#module load python/2.7
#rail-rna go local -x $IGENOME_HUMAN_BOWTIEINDEX $IGENOME_HUMAN_BOWTIE2INDEX -m $study_name/$sample_name/manifest -o $study_name/$sample_name/recount2 -d tsv,bw,jx -p 8 --scratch /mnt/work1/users/home2/zsaikhan/bhklab.home/temp_rail
#rm $study_name/$sample_name/*.sam
#date
