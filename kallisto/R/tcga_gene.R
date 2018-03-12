##create biobase expression object

#For aggregation of transcripts estimated counts/abundance to genes, it just adds them up to one number, no fancy calculation is included
#For length of the genes not sure how they are computed since it is not average length of transcripts as I checked
#So we can do it easily for UCSC/xena/Toil tcga and gtex
#keep in mind, in case of xena data, they are already transformed by log2(x+1)
##tcga
library(parallel)
options(stringsAsFactors=FALSE)
file.name <- "/mnt/work1/users/bhklab/Users/zhaleh/kallisto.gencode.v23/tcga_Kallisto_est_counts"
#tcga <- read.table(gzfile(file.name), header=T, row.names=1)
library(data.table)
message("Warning   [[[[file size: 14,813,269,044 alomist 15 GB]]]]")
#tcga <- fread(file.name) 
#rownames(tcga) <- tcga$sample   
incept <- ifelse(length(grep("counts", file.name)) > 0, 1, 0.001)
length(intersect(rownames(tcga), toil.transcripts$transcript_id))

tcga.gene.matrix <- matrix(0, ncol=ncol(tcga), nrow=length(unique(toil.transcripts$gene_id)))

colnames(tcga.gene.matrix) <- colnames(tcga)
rownames(tcga.gene.matrix) <- unique(toil.transcripts$gene_id)

ncore <- detectCores()
splitix <- parallel::splitIndices(nx=nrow(tcga.gene.matrix), ncl=ncore)
splitix <- splitix[sapply(splitix, length) > 0]

parallel::mclapply(splitix, function(x){
  geneid <- rownames(tcga.gene.matrix)[x]
  transids <- toil.transcripts[which(toil.transcripts$gene_id == geneid), "transcript_id"]
  apply(colnames(tcga.gene.matrix), function(sample){
    tcga.gene.matrix[geneid, sample] <- log2(sum(2^(as.numeric(tcga[transids, sample])) - incept) + incept)
  })
}, mc.cleanup=ncore)

for(geneid in rownames(tcga.gene.matrix)){
  transids <- toil.transcripts[which(toil.transcripts$gene_id == geneid), "transcript_id"]
  for(sample in colnames(tcga.gene.matrix)){
    tcga.gene.matrix[geneid, sample] <- log2(sum(2^(as.numeric(tcga[transids, sample])) - incept) + incept)
  }
}