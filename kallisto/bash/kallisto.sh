#!/bin/bash
#$ -S /bin/bash
#$ -N kallisto
#$ -o /mnt/work1/users/home2/zsaikhan/run/logs
#$ -e /mnt/work1/users/home2/zsaikhan/run/logs
module load hisat2/2.0.5
module load samtools/1.3.1
#module load stringtie/1.3.1c
##The 1.3.3b version of stringtie is installed in home directory and located in $HOME/bin
module load picard  
module load kallisto/0.43.1
module load igenome-human/hg38
## Usage
## CCLE: $1 study, $2 name of cell line, $3 bam file
## GNE: $1 study, $2 name of cell line, $3 input_fastq_1, $4 input_fastq_2
## UHN: $1 study, $2 name of sample, $3 input_fastq_1, $4 input_fastq_2
## GRAY: $1 study, $2 name of sample
# Navigate to file
cd  /mnt/work1/users/bhklab/Users/zhaleh/kallisto.hg38
if [ ! -d "$1/$2" ]; then
	mkdir -p $1/$2
fi

date
if [ "$1" == "PDX_DAVE" ]; then
    echo "Running kallisto"
    echo $1
    echo $2
    echo $3
    echo $4
   # hisat2 -p 12 --dta -x $HOME/Genome/GRCh38/Hisat/grch38_tran/grch38_tran -1 $3 -2 $4 -S $1/$2/Aligned.out.sam
   kallisto quant -t 8 -i /mnt/work1/users/bhklab/Users/zhaleh/hg38/kallisto_hg38.idx -o $1/$2 $3 $4
fi

if [ "$1" == "GNE" ]; then
    echo "Running kallisto"
    echo $1
    echo $2
    echo $3
    echo $4
   # hisat2 -p 12 --dta -x $HOME/Genome/GRCh38/Hisat/grch38_tran/grch38_tran -1 $3 -2 $4 -S $1/$2/Aligned.out.sam
   kallisto quant -t 8 -i /mnt/work1/users/bhklab/Users/zhaleh/hg38/kallisto_hg38.idx -o $1/$2 $3 $4
fi

if [ "$1" == "UHN" ]; then
    echo "Running kallisto"
    echo $1
    echo $2
    echo $3
    echo $4

   # hisat2 -p 12 --dta -x $HOME/Genome/GRCh38/Hisat/grch38_tran/grch38_tran -1 $3 -2 $4 -S $1/$2/Aligned.out.sam
   kallisto quant -t 8  -i /mnt/work1/users/bhklab/Users/zhaleh/hg38/kallisto_hg38.idx -o $1/$2 $3 $4

fi

if [ "$1" == "CCLE" ]; then
	echo "Running picard"
        echo $1
        echo $2
        echo $3

	java -Xmx16g -jar $picard_dir/picard.jar SamToFastq  I=$3 FASTQ=$1/$2/first.fastq SECOND_END_FASTQ=$1/$2/second.fastq
	echo "Running kallisto"
	#hisat2 -p 12 --dta -x $HOME/Genome/GRCh38/Hisat/grch38_tran/grch38_tran -1 $1/$2/first.fastq -2 $1/$2/second.fastq -S $1/$2/Aligned.out.sam
	kallisto quant -t 8  -i /mnt/work1/users/bhklab/Users/zhaleh/hg38/kallisto_hg38.idx -o $1/$2 $1/$2/first.fastq $1/$2/second.fastq
	rm $1/$2/first.fastq
	rm $1/$2/second.fastq
fi

if [ "$1" == "GRAY" ]; then
   	echo"Running fastq-dump"
        echo $1
        echo $2

	fastq-dump --split-files /mnt/work1/users/bhklab/Data/Gray/rnaseq/SRP026537/$2/$2.sra --outdir $1/$2/
	echo $2_1.fastq
	echo $2_2.fastq
	echo "Running kallisto"
	#hisat2 -p 12 --dta -x $HOME/Genome/GRCh38/Hisat/grch38_tran/grch38_tran -1 $1/$2/$2_1.fastq -2 $1/$2/$2_2.fastq -S $1/$2/Aligned.out.sam
	kallisto quant -t 8  -i /mnt/work1/users/bhklab/Users/zhaleh/hg38/kallisto_hg38.idx -o $1/$2 $1/$2/$2_1.fastq $1/$2/$2_2.fastq
	rm $1/$2/$2_1.fastq
	rm $1/$2/$2_2.fastq
fi
date

if [ "$1" == "GEO" ]; then
   	echo"Running fastq-dump"
        echo $1
        echo $2
        echo $3
        echo $4
    if [ ! -d "$1/$2/$3" ]; then
    	mkdir -p $1/$2/$3
    fi

	fastq-dump --split-files $4 --outdir $1/$2/$3
	echo $3_1.fastq
	echo $3_2.fastq
	echo "Running kallisto"
	#hisat2 -p 12 --dta -x $HOME/Genome/GRCh38/Hisat/grch38_tran/grch38_tran -1 $1/$2/$2_1.fastq -2 $1/$2/$2_2.fastq -S $1/$2/Aligned.out.sam
	kallisto quant -t 8  -i /mnt/work1/users/bhklab/Users/zhaleh/hg38/kallisto_hg38.idx -o $1/$2/$3 $1/$2/$3/$3_1.fastq $1/$2/$3/$3_2.fastq		
	rm $1/$2/$3/$3_1.fastq
	rm $1/$2/$3/$3_2.fastq
fi
date
#echo "Sort and Convert Sam file to bam"
#samtools sort -@ 8 -o $1/$2/Aligned.out.sorted.bam $1/$2/Aligned.out.sam

#date
#echo "Running stringtie"
#stringtie $1/$2/Aligned.out.sorted.bam -v -o $1/$2/stringtie_output.gtf -A $1/$2/gene_abund.tab -p 8 -G $HOME/Genome/GRCh38/Gencode/gencode.v26.annotation.gtf
#stringtie $1/$2/Aligned.out.sorted.bam -v -o $1/$2/ballgown/stringtie_output.gtf -e -B -p 8 -G $HOME/Genome/GRCh38/Gencode/gencode.v26.annotation.gtf
#stringtie $1/$2/Aligned.out.sorted.bam -v -o $1/$2/test/stringtie_output.gtf -p 8 -G $HOME/Genome/GRCh38/Gencode/gencode.v26.annotation.gtf
#date

#echo "Running Rail-rna"
#module load igenome-human/hg38
#module load python/2.7
#rail-rna go local -x $IGENOME_HUMAN_BOWTIEINDEX $IGENOME_HUMAN_BOWTIE2INDEX -m $1/$2/manifest -o $1/$2/recount2 -d tsv,bw,jx -p 8 --scratch /mnt/work1/users/home2/zsaikhan/bhklab.home/temp_rail
#rm $1/$2/*.sam
#date
