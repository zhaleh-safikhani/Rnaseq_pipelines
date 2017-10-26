options(stringsAsFactors=FALSE)
#install.packages("readr")

load("data/Gencode.v23.annotation.RData")
tx2gene <- as.data.frame(cbind("transcript"=toil.transcripts$transcript_id, "gene"=toil.transcripts$gene_id))

dir <- "PDX_DAVE"
files <- list.files(dir, recursive = TRUE, full.names = T)
resFiles <- grep("abundance.h5", files)
resFiles <- files[resFiles]
length(resFiles)
library(Biobase)
switch(dir,
       "CCLE"={
         names(resFiles) <- sapply(strsplit(resFiles, "\\/"), function(x){x[2]})
         library(PharmacoGx)
         load("CCLE_hs.RData")
         pp <- pData(CCLE@molecularProfiles$rnaseq)
       },
       "GNE"={
         names(resFiles) <- sapply(strsplit(resFiles, "\\/"), function(x){x[2]})
         library(PharmacoGx)
         load("gCSI_hs.RData")
         pp <- pData(gCSI@molecularProfiles$rnaseq)
       },
       "UHN"={
         names(resFiles) <- sapply(strsplit(resFiles, "\\/"), function(x){x[2]})
         library(PharmacoGx)
           load("UHN_hs.RData", verbose=T)
           pp <- pData(UHN@molecularProfiles$rnaseq)
           pp[which(pp$cellid == "AU655"), "cellid"] <- "AU565"
           UHN@cell <- UHN@cell[-which(rownames(UHN@cell)== "AU655"),]
           UHN@curation$cell <- UHN@cell[-which(rownames(UHN@curation$cell)== "AU655"),]
           UHN@curation$tissue <- UHN@cell[-which(rownames(UHN@curation$tissue)== "AU655"),]
           load("UHNBreast.RData", verbose=T)
           UHNBreast@cell <- UHN@cell
           UHN <- UHNBreast
           UHN@sensitivity$n <- rbind(UHN@sensitivity$n, "SW 527"=rep(0,ncol(UHN@sensitivity$n))) 
       },
       "GRAY"={
         names(resFiles) <- sapply(strsplit(resFiles, "\\/"), function(x){x[2]})
         library(PharmacoGx)
         load("GRAY_hs.RData")
         pp <- pData(GRAY@molecularProfiles$rnaseq)
       },
       "GEO/GSE79668"={
         runs <- sapply(strsplit(resFiles, "\\/"), function(x){x[3]})
         #source("http://bioconductor.org/biocLite.R")
         #biocLite("GEOquery")
         library(GEOquery)
         gg <- getGEO("GSE79668")
         pp <- pData(gg$GSE79668_series_matrix.txt.gz)
         
         downloader::download(url = "http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=SRP072492", destfile="SRP072492.tab")
         ss <- read.csv(file ="SRP072492.tab", header=T)
         
         names(resFiles) <- ss[match(ss$Run, runs), "SampleName"]
         pp$run_id <- ss[match(rownames(pp), samples), "Run"]
         #SRR to GSM mapping solution2
         #biocLite("SRAdb")
         #library(SRAdb)
         #sqlfile <- getSRAdbFile()
         #sra_con <- dbConnect(SQLite(),sqlfile)
         #res <- dbGetQuery(sra_con, "select * from sra_ft where run_accession='SRR097786'")
         
         #SRR to GSM mapping solution3
         #downloader::download(url = "ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab", destfile="~/Google Drive/CCLEDrug/data/SRA_Accessions.tab")
         #SRA.accession <- read.table("~/Google Drive/CCLEDrug/data/SRA_Accessions.tab")
       },
       "PDX_DAVE"={
         names(resFiles) <- sapply(strsplit(resFiles, "\\/"), function(x){x[2]})
         pp <- as.data.frame(matrix(NA, ncol=1, nrow=length(resFiles)))
         rownames(pp) <- pp[,1] <- names(resFiles)
         colnames(pp) <- "sample_id"
       })


library(readr)
library(tximport)
txi <- tximport(resFiles, type="kallisto", tx2gene=tx2gene)
head(txi$counts[,1:5])
dim(txi$counts)

xx <- txi$abundance
gene.exp <- Biobase::ExpressionSet(log2(xx + 0.001))
fData(gene.exp) <- toil.genes[featureNames(gene.exp),]
pData(gene.exp) <- pp[sampleNames(gene.exp),]
annotation(gene.exp) <- "rnaseq"

xx <- txi$counts
gene.count <- Biobase::ExpressionSet(log2(xx + 1))
fData(gene.count) <- toil.genes[featureNames(gene.count),]
pData(gene.count) <- pp[sampleNames(gene.count),]
annotation(gene.count) <- "rnaseq"

txii <- tximport(resFiles, type="kallisto", txOut=T)

xx <- txii$abundance
transcript.exp <- Biobase::ExpressionSet(log2(xx[,1:length(resFiles)] + 0.001))
fData(transcript.exp) <- toil.transcripts[featureNames(transcript.exp),]
pData(transcript.exp) <- pp[sampleNames(transcript.exp),]
annotation(transcript.exp) <- "isoforms"

xx <- txii$counts
transcript.count <- Biobase::ExpressionSet(log2(xx[,1:length(resFiles)] + 1))
fData(transcript.count) <- toil.transcripts[featureNames(transcript.count),]
pData(transcript.count) <- pp[sampleNames(transcript.count),]
annotation(transcript.count) <- "isoforms"

switch(dir,
       "CCLE"={
         CCLE@molecularProfiles$rnaseq <- gene.exp
         CCLE@molecularProfiles$rnaseq.counts <- gene.count
         CCLE@molecularProfiles$isoforms <- transcript.exp
         CCLE@molecularProfiles$isoforms.counts <- transcript.count
         save(CCLE, file="CCLE_kallisto.RData")
         ccle.rnaseq.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=CCLE, mDataType="rnaseq", fill.missing=FALSE)))
         
      },
       "GNE"={
         gCSI@molecularProfiles$rnaseq <- gene.exp
         gCSI@molecularProfiles$rnaseq.counts <- gene.count
         gCSI@molecularProfiles$isoforms <- transcript.exp
         gCSI@molecularProfiles$isoforms.counts <- transcript.count
         save(gCSI, file="gCSI_kallisto.RData")
         gcsi.rnaseq.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=gCSI, mDataType="rnaseq", fill.missing=FALSE)))
         
       },
       "UHN"={
         UHN@molecularProfiles$rnaseq <- gene.exp
         UHN@molecularProfiles$rnaseq.counts <- gene.count
         UHN@molecularProfiles$isoforms <- transcript.exp
         UHN@molecularProfiles$isoforms.counts <- transcript.count
         save(UHN, file="UHN_kallisto.RData")
         uhn.rnaseq.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=UHN, mDataType="rnaseq", fill.missing=FALSE)))
       },
       "GRAY"={
         GRAY@molecularProfiles$rnaseq <- gene.exp
         GRAY@molecularProfiles$rnaseq.counts <- gene.count
         GRAY@molecularProfiles$isoforms <- transcript.exp
         GRAY@molecularProfiles$isoforms.counts <- transcript.count
         save(GRAY, file="GRAY_kallisto.RData")
       },
       "GEO/GSE79668"={
         save(gene.exp, gene.count, transcript.exp, transcript.count, file="GSE79668.kallisto.hg38.RData")
       },
      "PDX_DAVE"={
        save(gene.exp, gene.count, transcript.exp, transcript.count, file="gene.exp.kallisto.hg38.RData")
        })



