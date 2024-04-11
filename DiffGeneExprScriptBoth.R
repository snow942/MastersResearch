###Import libraries###
library (GenomicFeatures) 
library (edgeR)
library (limma)
library (gplots)
library (ggplot2)
library (tximport)
library (readr)
library (jsonlite)
library (matrixStats)
library (cowplot)
library (matrixStats)
library (DT)
library (plotly)
library (gt)
library (RColorBrewer)
library (heatmaply)

setwd("C:/Users/snow4/Desktop/kallisto")

#Kallisto previously run via command line


###Making annotation file###
#Protocol provided by Dr. Geneva, will post directly below
#


###Kallisto data import###
dir <- "./"
list.files (file.path (dir, "kallisto"))
samples <- read.table (file.path (dir, "samplesBoth.txt"), header=TRUE)
files <- file.path (dir, "kallisto", 
                    samples$sample, 
                    "abundance.tsv")
names (files) <- paste0 ("sample", 1:5)


###Making txdb and tx2gene of non-model from Ensembl gtf file###
txdb <- makeTxDbFromGFF ("Adineta_vaga.AMS_PRJEB1171_v1.50.gtf", format=c("gtf"))
k <- keys (txdb, keytype= "TXNAME")
tx2gene <- select (txdb, 
                   k, 
                   columns = "GENEID", 
                   keytype = "TXNAME")

library (tidyverse)
library (dplyr)
library (rhdf5)
#Packages need to be imported after to ensure select () functions properly

txi <- tximport (files, 
                 type="kallisto", 
                 tx2gene=tx2gene, 
                 txOut = FALSE, 
                 countsFromAbundance = "lengthScaledTPM", 
                 ignoreTxVersion = TRUE)


###Normalization and filtering of data###
TPM <- txi$abundance
Counts <- txi$counts
sampleLabels <- samples$sample
DGEList <- DGEList (Counts)
CPM <- cpm (DGEList)

log2.cpm <- cpm (DGEList, log = TRUE)
log2.cpm.df <- as_tibble (log2.cpm, rownames = "GeneName")
colnames (log2.cpm.df) <- c("GeneName", sampleLabels)
log2.cpm.df.pivot <- pivot_longer (log2.cpm.df,
                                   cols = Lane1Lang:SRR7962071,
                                   names_to = "Samples",
                                   values_to = "Expression")

ggplot (log2.cpm.df.pivot) +
  aes (x = Samples, y = Expression, fill = Samples) +
  geom_violin (trim = FALSE, show.legend = FALSE) +
  labs (y = "Log2 Expression", x = "Sample",
        title = "Langjokull and Vatnajokull Phiodina sp. Log2 Counts Per Million (CPM)",
        subtitle = "Unfiltered, non-normalized",
        caption = "Ryan Snow") + 
  theme_classic()


Keep <- rowSums (CPM > 1) > 4
#Filter to keep genes with CPM that is greater than 1 in all of the samples
DGEList.Filtered <- DGEList [Keep,]

log2.cpm.filtered <- cpm (DGEList.Filtered, log = TRUE)
log2.cpm.filtered.df <- as_tibble (log2.cpm.filtered, rownames = "GeneName")
colnames (log2.cpm.filtered.df) <- c("GeneName", sampleLabels)
log2.cpm.filtered.df.pivot <- pivot_longer (log2.cpm.filtered.df,
                                   cols = Lane1Lang:SRR7962071,
                                   names_to = "Samples",
                                   values_to = "Expression")

ggplot (log2.cpm.filtered.df.pivot) +
  aes (x = Samples, y = Expression, fill = Samples) +
  geom_violin (trim = FALSE, show.legend = FALSE) +
  stat_summary (fun = "median",
                geom = "point",
                shape = 95,
                size  = 10,
                color = "black",
                show.legend = FALSE) +
  labs (y = "Log2 Expression", x = "Sample",
        title = "Langjokull and Vatnajokull Philodina sp. Log2 Counts Per Million (CPM)",
        subtitle = "Filtered, non-normalized",
        caption = "Ryan Snow") + 
  theme_classic()

DGEList.Filtered.norm <- calcNormFactors (DGEList.Filtered, method = "TMM")
#TMM normalization of data

log2.cpm.filtered.norm <- cpm (DGEList.Filtered.norm, log = TRUE)
log2.cpm.filtered.norm.df <- as_tibble (log2.cpm.filtered.norm, rownames = "GeneName")
colnames (log2.cpm.filtered.norm.df) <- c("GeneName", sampleLabels)
log2.cpm.filtered.norm.df.pivot <- pivot_longer (log2.cpm.filtered.norm.df,
                                            cols = Lane1Lang:SRR7962071,
                                            names_to = "Samples",
                                            values_to = "Expression")

ggplot (log2.cpm.filtered.norm.df.pivot) +
  aes (x = Samples, y = Expression, fill = Samples) +
  geom_violin (trim = FALSE, show.legend = FALSE) +
  stat_summary (fun = "median",
                geom = "point",
                shape = 95,
                size  = 10,
                color = "black",
                show.legend = FALSE) +
  labs (y = "Log2 Expression", x = "Sample",
        title = "Langjokull and Vatnajokull Philodina sp. Log2 Counts Per Million (CPM)",
        subtitle = "Filtered, TMM normalized",
        caption = "Ryan Snow") + 
  theme_classic()

###PCA###
group <- factor (samples$condition)
distance <- dist (t (log2.cpm.filtered.norm), method = "maximum")
clusters <- hclust (distance, method = "complete")
plot (clusters, labels = sampleLabels)
pca.res <- prcomp (t (log2.cpm.filtered.norm), 
                   scale. = F, 
                   retx = T)

screeplot (pca.res)
#summary (pca.res)
pc.var <- pca.res$sdev^2
pc.per <- round (pc.var/sum (pc.var)*100, 1)
pca.res.df <- as_tibble (pca.res$x)

ggplot (pca.res.df) +
  aes(x = PC1, y = PC2, label = sampleLabels) +
  geom_point (size = 2) +
  geom_label() +
  #stat_ellipse () +
  xlab (paste0 ("PC1 (", pc.per[1], "%)")) +
  ylab (paste0 ("PC2 (", pc.per[2], "%)")) +
  labs (title = "PCA Plot",
        caption = "Ryan Snow") +
#  coord_fixed () +
  theme_classic()

pca.res.df <- pca.res$x[,1:4] %>%
  as_tibble() %>%
  add_column (sample = sampleLabels,
              group = group)
pca.pivot <- pivot_longer (pca.res.df,
                           cols = PC1:PC4,
                           names_to = "PC",
                           values_to = "Loadings")

ggplot (pca.pivot) +
  aes (x = sample, y = Loadings, fill = group) +
  geom_bar (stat = "Identity") +
  facet_wrap (~PC) +
  labs (title = "PCA 'Small Multiples' Plot",
        caption = "Ryan Snow") +
  theme_classic() +
  coord_flip()

data.df <- log2.cpm.filtered.norm.df %>%
  mutate (Psychro.AVG = (Lane1Lang + Lane2Lang + Lane1Vat + Lane2Vat)/4,
          Meso.AVG = (SRR7962071 + SRR7962070 + SRR7962068)/3,
          LogFC = (Psychro.AVG - Meso.AVG)) %>%
  mutate_if (is.numeric, round, 2)

data.sort <- data.df %>%
  dplyr::arrange (desc (LogFC)) %>%
  dplyr::select (GeneName, LogFC)


###Differential Expression Analysis###
design <- model.matrix (~0 + group)
colnames (design) <- levels (group)

v.DEGList.filtered.norm <- voom (DGEList.Filtered.norm,
                                 design,
                                 plot = TRUE)
fit <- lmFit (v.DEGList.filtered.norm, design)
contrast.matrix <- makeContrasts (ColdAdaptation = Psychro - Meso,
                                  levels = design)
fits <- contrasts.fit (fit, contrast.matrix)
ebFit <- eBayes (fits)
TopHits <- topTable (ebFit, 
                        adjust = "BH", 
                        coef = 1, 
                        number = 50000,
                        sort.by = "logFC")

TopHits.df <- TopHits %>%
  as_tibble (rownames = "GeneName") 
gt (TopHits.df)

VolcanoPlot <- ggplot (TopHits.df) +
  aes (y = -log10(adj.P.Val), x = logFC, text = paste("Symbol:", GeneName)) +
  geom_point (size = 2) +
  labs (title = "Langjokull and Vatnajokull Philodina sp.",
        subtitle = "Differentially expressed targets",
        caption = "Ryan Snow") +
  theme_classic()

VolcanoPlot
ggplotly (VolcanoPlot)

results <- decideTests (ebFit, method = "global", adjust.method = "BH", p.value = 0.01, lfc = 2)
summary (results)

vennDiagram (results, include = "both")

head(v.DEGList.filtered.norm$E)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels

diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
diffGenesIncrease <- v.DEGList.filtered.norm$E[results[,1] ==1,]
diffGenesDecrease <- v.DEGList.filtered.norm$E[results[,1] ==-1,]

diffGenes.df <- as_tibble(diffGenes, rownames = "GeneName")
diffGenesIncrease.df <- as_tibble(diffGenesIncrease, rownames = "GeneName")
diffGenesDecrease.df <- as_tibble(diffGenesDecrease, rownames = "GeneName")
datatable(diffGenes.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Table 1: Genes differentially expressed',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("25", "50", "100", "250"))) %>%
  formatRound(columns=c(2:6), digits=2)

write_tsv (diffGenes.df,"DiffGenesBoth.txt")
write_tsv (diffGenesIncrease.df,"DiffGenesIncreaseBoth.txt")
write_tsv (diffGenesDecrease.df,"DiffGenesDecreaseBoth.txt")

NewTopHits <- topTable (ebFit, adjust = "BH", coef = 1, number = 100, sort.by = "logFC", p.value = 0.01)
NewTopHits.df <- as_tibble (NewTopHits, rownames = "GeneName")
write_tsv (NewTopHits.df,"NewTopHitsBoth.txt")

###Annotation###
#DiffGenes.txt was merged with aforementioned annotation file in Excel
#Merged by loading files as data queries, matching via "GeneName" (i.e. "GSADVG..."), then using the merge function provided in Excel


###Functional annotation###
#DiffGenesIncrease.txt, DiffGenesDecrease.txt, and NewTopHits.txt submitted to gProfiler web server for GO annotations
  #Manhattan plot(s) created
    #Default options

###Heatmap generation###
Color <- brewer.pal (name = "PiYG", n=11)
clustRows <- hclust (as.dist (1-cor (t (diffGenes), 
                                     method = "pearson")),
                     method = "complete")
clustColumns <- hclust (as.dist (1-cor (diffGenes, 
                                        method = "pearson")),
                        method = "complete")
module.assign <- cutree (clustRows, k = 2)
module.color <- rainbow (length (unique (module.assign)), start = 0.1, end = 0.9)
module.color <- module.color[as.vector (module.assign)]

heatmap.2 (diffGenes,
           Rowv = as.dendrogram (clustRows),
           Colv = as.dendrogram (clustColumns),
           RowSideColors = module.color,
           col = Color,
           scale = "row",
           labRow = NA,
           density.info = "none",
           trace = "none",
           cexRow = 1,
           cexCol = 1,
           margins = c(8, 20))

