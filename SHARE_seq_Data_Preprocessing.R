library(Signac)
library(Seurat)
library(ggplot2)
library(EnsDb.Mmusculus.v79)
library(magrittr)

#### Load Data####
#---> load processed data matrices for each assay----
pathway <- ""
rna <- Read10X(paste(pathway, "SHARE_RNA-seq", sep = ""), gene.column = 1)
atac <- Read10X(paste(pathway, "SHARE_ATAC-seq", sep = ""), gene.column = 1)

# create a Seurat object and add the assays
share <- CreateSeuratObject(counts = rna)

share[['ATAC']] <- CreateChromatinAssay(
  counts = atac,
  sep = c(":", "-"),
)

remove(rna, atac)
gc()

#----> Annotation ----
# extract gene annotations from EnsDb
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
# change to UCSC style
genome(annotation) <- "mm10"
annotation <- renameSeqlevels(annotation, mapSeqlevels(seqlevels(annotation), "UCSC"))
#seqlevelsStyle(annotation) <- "UCSC"
Annotation(share[["ATAC"]]) <- annotation

#----> Celltype annotation ----
celltype_df <- read.table(paste(pathway, "GSM4156597_skin_celltype.txt", sep = ""), header = TRUE, sep = "\t")
colOrd <- rownames(FetchData(share,"ident"))
celltype_df <- celltype_df[match(noquote(colOrd), celltype_df$atac.bc),]

share@meta.data <- cbind(share@meta.data, celltype_df$celltype)
colnames(share@meta.data)[which(names(share@meta.data) == "celltype_df$celltype")] <- "celltype"

# Cell type summary df
celltype.summary <- data.frame(summary(factor(share@meta.data[[6]])))
colnames(celltype.summary)[1] <- "counts"
celltype.summary <- cbind(Celltype = rownames(celltype.summary), celltype.summary)
celltype.summary <- celltype.summary[order(-celltype.summary$counts),]
rownames(celltype.summary) <- 1:nrow(celltype.summary)

# Filter
#----> RNA QC ----
# Remove low quality cells
DefaultAssay(share) <- "RNA"
share <- subset(
  x = share,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 2500
)
# Remove genes on sex chromosome
SexChro_genes <- read.table(file = paste(pathway, "sex_chromo.txt", sep = ""), header = FALSE, sep = "\t")
colnames(SexChro_genes) <- c("gene")

DefaultAssay(share) <- "RNA"
counts <- GetAssayData(share, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% SexChro_genes$gene)),]
share@assays$RNA <- subset(x = share@assays$RNA, features = rownames(counts))

#----> ATAC QC ----
# Remove peaks on sex chromosome
peak_region <- data.frame(share@assays$ATAC@ranges)
counts <- GetAssayData(share, assay = "ATAC")
counts <- counts[-which(peak_region$seqnames == c("chrX")),]

# Binarize data
counts <- BinarizeCounts(counts)

# Remove too low & high peak expression
peaks_rowsum <- rowSums(counts)
peaks_toRemove1 <- peaks_rowsum[which(peaks_rowsum < 5)]
peaks_toRemove2 <- peaks_rowsum[which(peaks_rowsum > dim(counts)[2] * 0.1)]
peaks_toRemove <- c(names(peaks_toRemove1), names(peaks_toRemove2))
counts <- counts[-(which(rownames(counts) %in% peaks_toRemove)),]

share@assays$ATAC <- subset(x = share@assays$ATAC, features = rownames(counts))

#### Normalization ####
DefaultAssay(share) <- "RNA"
share <- NormalizeData(share, normalization.method = "LogNormalize", scale.factor = 1)
share <- ScaleData(share)

#### Convert Seurat to Scanpy ####
#remotes::install_github("pmbio/MuDataSeurat")
library(MuDataSeurat)
MuDataSeurat::WriteH5MU(share, paste(pathway, "share.h5mu", sep = ""))