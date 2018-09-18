#!/bin/bash
#$ -S /bin/bash
#$ -N STAR_genomeGenerate
#$ -o $HOME/run/logs
#$ -e $HOME/run/logs

# Load modules
module load STAR/2.4.2a
module load igenome-human/GRCh38
module load picard/1.9.1

# Hard-coded variables
GENOMEDIR=$1 # the output direcotory where the indices are located
FIRST_FASTQ=$2
SECOND_FASTQ=$3
OUTPUT=$4 # read length - 1

mkdir -p $OUTPUT'/tmp'

# --readFilesCommand if fastq files are gzipped: readFilesCommand gunzip -c
# --genomeSAsparseD default value is 1 use bigger numbers than 1 to reduce ram consumption
# --twopassMode Basic To run STAR 2-pass mapping for each sample separately
# --outSAMprimaryFlag AllBestScore all alignments with the best score are primary; OneBestScore only one alignment with the best score is primary
# --outFilterIntronMotifs RemoveNoncanonical filter out alignments that contain non-canonical junctions
# --outSAMtype BAM SortedByCoordinate output sorted by coordinate Aligned.sortedByCoord.out.bam file, similar to samtools sort command.
# --chimSegmentMin To switch on detection of chimeric (fusion) alignments (in addition to normal mapping), 
###--chimSegmentMin should be set to a positive value
###For example, if you have 2x75 reads and used --chimSegmentMin 20, a chimeric alignment with 130b on one chromosome and 20b on the other will be output, 
###while 135 + 15 won’t be.
#--chimOutType SeparateSAMold output old SAM into separate Chimeric.out.sam file
#--quantMode TranscriptomeSAM output SAM/BAM alignments to transcriptome into a separate file; With --quantMode TranscriptomeSAM GeneCounts, and get both the Aligned.toTranscriptome.out.bam
#and ReadsPerGene.out.tab outputs.
#--outWigNorm default: RPM string: type of normalization for the signal; RPM reads per million of mapped reads; None no normalization, ”raw” counts
#--limitIObufferSize default: 150000000 int>0: max available buffers size (bytes) for input/output, per thread
#--limitBAMsortRAM default: 0 int>=0: maximum available RAM for sorting BAM. . If =0, it will be set to the genome index size. 0 value can only be used with –genomeLoad NoSharedMemory option
# --outReadsUnmapped Fastx output in separate fasta/fastq files, Unmapped.out.mate1/2
# --outSAMunmapped Within output unmapped reads within the main SAM file (i.e. Aligned.out.sam)
# --outBAMsortingThreadN default: 0 int: >=0: number of threads for BAM sorting. 0 will default to min(6,–runThreadN).

date
echo "Running star"

STAR --runMode alignReads \
--genomeDir $GENOMEDIR \
--readFilesIn $FIRST_FASTQ $SECOND_FASTQ \
--readFilesCommand gunzip -c \
--runThreadN 6 \
--genomeSAsparseD 2 \
--twopassMode Basic \
--outFileNamePrefix $OUTPUT \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--quantMode TranscriptomeSAM \
--outSAMattributes NH HI AS NM MD \
--outFilterType BySJout \
--outFilterMultimapNmax 20 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.04 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--sjdbScore 1
date
echo "Star Complete"
