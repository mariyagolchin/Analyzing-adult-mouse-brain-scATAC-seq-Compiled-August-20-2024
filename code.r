
# https://stuartlab.org/signac/articles/mouse_brain_vignette

setwd("/home/golchinpour/projects/stuartlab_projs/1-scATAC-seq_MouseBrain_DNA_Motif_Analysis")
library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
################### Pre-processing workflow ###################
counts <- Read10X_h5("Mouse Brain data/atac_v1_adult_brain_fresh_5k_filtered_peak_bc_matrix.h5")
# > head(counts[1:5,1:3])
# 5 x 3 sparse Matrix of class "dgCMatrix"
#                      AAACGAAAGTAATCAG-1 AAACGAACACGCTGTG-1 AAACGAATCCTGGGAC-1
# chr1:3094708-3095552                  .                  .                  .
# chr1:3119414-3121782                  2                  .                  .
# chr1:3204809-3205178                  .                  .                  .
# chr1:3217330-3217359                  .                  .                  .
# chr1:3228123-3228463                  2                  .                  .
# > 

metadata <- read.csv(
  file = "Mouse Brain data/atac_v1_adult_brain_fresh_5k_singlecell.csv",
  header = TRUE,
  row.names = 1


brain_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = './Mouse Brain data/atac_v1_adult_brain_fresh_5k_fragments.tsv.gz',
  min.cells = 1
)

brain <- CreateSeuratObject(
  counts = brain_assay,
  assay = 'peaks',
  project = 'ATAC',
  meta.data = metadata
)

# add gene annotations to the brain object for the mouse genome.
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"

# add the gene information to the object
Annotation(brain) <- annotations

Annotation(brain)
> brain$
brain$orig.ident                        brain$lowmapq                           brain$enhancer_region_fragments
brain$nCount_peaks                      brain$mitochondrial                     brain$promoter_region_fragments
brain$nFeature_peaks                    brain$passed_filters                    brain$on_target_fragments
brain$total                             brain$cell_id                           brain$blacklist_region_fragments
brain$duplicate                         brain$is__cell_barcode                  brain$peak_region_fragments
brain$chimeric                          brain$TSS_fragments                     brain$peak_region_cutsites
brain$unmapped                          brain$DNase_sensitive_region_fragments  
> brain$

> brain@
brain@assays        brain@active.ident  brain@reductions    brain@misc          brain@tools         
brain@meta.data     brain@graphs        brain@images        brain@version       
brain@active.assay  brain@neighbors     brain@project.name  brain@commands      
> brain@

# Computing QC Metrics
# Next we compute some useful per-cell QC metrics.

brain <- NucleosomeSignal(object = brain)

brain$nucleosome_group <- ifelse(brain$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
fraghist<- FragmentHistogram(object = brain, group.by = 'nucleosome_group', region = 'chr1-1-10000000')

library(ggplot2)

# Save the plot to a file
ggsave(
  filename = "2-plots/1-FragmentHistogram.png",  # Specify the output file name
  plot = fraghist,           # The plot object to save
  device = "png",                 # File format; change to "pdf" for PDF
  width = 10,                     # Width of the plot in inches
  height = 6,                     # Height of the plot in inches
  dpi = 300                       # Resolution in dots per inch
)

brain <- TSSEnrichment(brain) # use index fragment file
brain$pct_reads_in_peaks <- brain$peak_region_fragments / brain$passed_filters * 100
brain$blacklist_ratio <- brain$blacklist_region_fragments / brain$peak_region_fragments

vlnplot<- VlnPlot(
  object = brain,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
ggsave(
  filename = "2-plots/2-vlnplot.png",  # Specify the output file name
  plot = vlnplot,           # The plot object to save
  device = "png",                 # File format; change to "pdf" for PDF
  width = 10,                     # Width of the plot in inches
  height = 6,                     # Height of the plot in inches
  dpi = 300                       # Resolution in dots per inch
)

brain <- subset(
  x = brain,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.025 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

vlnplot_rmoutliers <- VlnPlot(
  object = brain,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
ggsave(
  filename = "2-plots/2-vlnplot_vlnplot_rmoutliers.png",  # Specify the output file name
  plot = vlnplot_rmoutliers,           # The plot object to save
  device = "png",                 # File format; change to "pdf" for PDF
  width = 10,                     # Width of the plot in inches
  height = 6,                     # Height of the plot in inches
  dpi = 300                       # Resolution in dots per inch
)
# Normalization and linear dimensional reduction
brain <- RunTFIDF(brain)
> brain <- FindTopFeatures(brain, min.cutoff = 'q0')
> brain <- RunSVD(object = brain)
depcor<- DepthCor(brain)

ggsave(
  filename = "2-plots/3-depcor.png",  # Specify the output file name
  plot = depcor,           # The plot object to save
  device = "png",                 # File format; change to "pdf" for PDF
  width = 10,                     # Width of the plot in inches
  height = 6,                     # Height of the plot in inches
  dpi = 300                       # Resolution in dots per inch
)

brain <- RunUMAP(
  object = brain,
  reduction = 'lsi',
  dims = 2:30
)
brain <- FindNeighbors(
  object = brain,
  reduction = 'lsi',
  dims = 2:30
)
brain <- FindClusters(
  object = brain,
  algorithm = 3,
  resolution = 1.2,
  verbose = FALSE
)

dimplot<- DimPlot(object = brain, label = TRUE) 

ggsave(
  filename = "2-plots/4-dimplot.png",  # Specify the output file name
  plot = dimplot,           # The plot object to save
  device = "png",                 # File format; change to "pdf" for PDF
  width = 10,                     # Width of the plot in inches
  height = 6,                     # Height of the plot in inches
  dpi = 300                       # Resolution in dots per inch
)

# ===================================================================
#                 Try to label clusters
# ===================================================================

# Create a gene activity matrix
# compute gene activities
gene.activities <- GeneActivity(brain) # use index fragment file
> head(gene.activities[1:6,1:3])
# 6 x 3 sparse Matrix of class "dgCMatrix"
#               AAACGAAAGTAATCAG-1 AAACGAACACGCTGTG-1 AAACGAATCCTGGGAC-1
# Hnf4g                          .                  .                  .
# Zfhx4                          2                  .                  .
# Pex2                           1                  1                  1
# UBC                            .                  .                  .
# 1700008P02Rik                  .                  .                  .
# Pkia                           4                  3                  2
# > 

brain[['RNA']] <- CreateAssayObject(counts = gene.activities)
brain <- NormalizeData(
  object = brain,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(brain$nCount_RNA)
)

DefaultAssay(brain) <- 'RNA'

fplot<- FeaturePlot(
  object = brain,
  features = c('Sst','Pvalb',"Gad2","Neurod6","Rorb","Syt6"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
ggsave(
  filename = "2-plots/5-fplot.png",  # Specify the output file name
  plot = fplot,           # The plot object to save
  device = "png",                 # File format; change to "pdf" for PDF
  width = 10,                     # Width of the plot in inches
  height = 6,                     # Height of the plot in inches
  dpi = 300                       # Resolution in dots per inch
)

####################################################
          Integrating with scRNA-seq data
####################################################
# Load the pre-processed scRNA-seq data
allen_rna <- readRDS("Mouse Brain data/allen_brain.rds")
allen_rna <- UpdateSeuratObject(allen_rna)

allen_rna <- FindVariableFeatures(
  object = allen_rna,
  nfeatures = 5000
)

transfer.anchors <- FindTransferAnchors(
  reference = allen_rna,
  query = brain,
  reduction = 'cca',
  dims = 1:30
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = allen_rna$subclass,
  weight.reduction = brain[['lsi']],
  dims = 2:30
)

brain <- AddMetaData(object = brain, metadata = predicted.labels)

plot1 <- DimPlot(allen_rna, group.by = 'subclass', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
plot2 <- DimPlot(brain, group.by = 'predicted.id', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
combineplot<- plot1 + plot2

ggsave(
  filename = "2-plots/6-combineplot.png",  # Specify the output file name
  plot = combineplot,           # The plot object to save
  device = "png",                 # File format; change to "pdf" for PDF
  width = 10,                     # Width of the plot in inches
  height = 6,                     # Height of the plot in inches
  dpi = 300                       # Resolution in dots per inch
)

# Find differentially accessible peaks between clusters
DefaultAssay(brain) <- 'peaks'
Idents(brain) <- "predicted.id"

da_peaks <- FindMarkers(
  object = brain,
  ident.1 = c("L2/3 IT"), 
  ident.2 = c("L4", "L5 IT", "L6 IT"),
  test.use = 'wilcox',
  min.pct = 0.1
)

head(da_peaks)

plot1 <- VlnPlot(
  object = brain,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  idents = c("L4","L5 IT","L2/3 IT")
)
plot2 <- FeaturePlot(
  object = brain,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  max.cutoff = 'q95'
)
combineplot<- plot1 + plot2

ggsave(
  filename = "2-plots/7-combineplot2.png",  # Specify the output file name
  plot = combineplot,           # The plot object to save
  device = "png",                 # File format; change to "pdf" for PDF
  width = 10,                     # Width of the plot in inches
  height = 6,                     # Height of the plot in inches
  dpi = 300                       # Resolution in dots per inch
)

open_l23 <- rownames(da_peaks[da_peaks$avg_log2FC > 3, ])
open_l456 <- rownames(da_peaks[da_peaks$avg_log2FC < 3, ])
closest_l23 <- ClosestFeature(brain, open_l23)
closest_l456 <- ClosestFeature(brain, open_l456)
head(closest_l23)
head(closest_l456)

# Plotting genomic regions

(which(table(Idents(brain)) > 50))
#      L6 IT    L2/3 IT Macrophage      Oligo      Pvalb      Lamp5      L5 IT 
#          1          2          3          4          5          6          7 
#        Vip      Astro         L4        Sst      L6 CT        L6b      L5 PT 
#          8          9         10         11         12         13         14 
# > 
idents.plot <- names(which(table(Idents(brain)) > 50))
> idents.plot
 [1] "L6 IT"      "L2/3 IT"    "Macrophage" "Oligo"      "Pvalb"     
 [6] "Lamp5"      "L5 IT"      "Vip"        "Astro"      "L4"        
[11] "Sst"        "L6 CT"      "L6b"        "L5 PT"     
> 
brain <- SortIdents(brain)
covplot<- CoveragePlot(
  object = brain,
  region = c("Neurod6", "Gad2"),
  idents = idents.plot,
  extend.upstream = 1000,
  extend.downstream = 1000,
  ncol = 1
)

ggsave(
  filename = "2-plots/8-covplot.png",  # Specify the output file name
  plot = covplot,           # The plot object to save
  device = "png",                 # File format; change to "pdf" for PDF
  width = 10,                     # Width of the plot in inches
  height = 6,                     # Height of the plot in inches
  dpi = 300                       # Resolution in dots per inch
)