#!/bin/bash
#$ -S /bin/bash
#$ -N STAR_genomeGenerate
#$ -o $HOME/run/logs
#$ -e $HOME/run/logs

module load igenome-human/GRCh38
module load STAR/2.4.0j

GENOMEDIR=$1 # the output direcotory where the indices are located
ANNOTATION_FILE=$2 # annotation file
READ_LENGHT=$3 # read length - 1

# sjdbOverhang = read lenght -1
# genomeChrBinNbits = min(18, log2(GenomeLength/NumberOfReferences) to reduce ram consumption
# genomeSAsparseD default value is 1 use bigger numbers than 1 to reduce ram consumption
# limitGenomeGenerateRAM default: 31000000000 maximum available RAM (bytes) for genome generation

date
echo "Running star in genome generate mode to index the annotation"
# Run Genome Generate
STAR --runMode genomeGenerate \
--runThreadN 12 \
--genomeDir $GENOMEDIR \
--genomeFastaFiles $REF \
--sjdbGTFfile $ANNOTATION_FILE \
--sjdbOverhang $READ_LENGHT \
--genomeChrBinNbits 18 \
--genomeSAsparseD 2 \
--limitGenomeGenerateRAM 40000000000 \
--outTmpDir $GENOMEDIR/GenomeTmp


