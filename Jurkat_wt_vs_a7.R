###Jurkat wt vs a7
###load packages ----
library(tidyverse)
library(tximport)
library(ensembldb) 
library(EnsDb.Hsapiens.v86)
library(edgeR)
library(matrixStats)
library(cowplot)
library(DT)
library(gt)
library(plotly)
library(limma)
library(RColorBrewer)
library(GSEABase) 
library(Biobase) 
library(GSVA) 
library(gprofiler2) 
library(clusterProfiler) 
library(msigdbr) 
library(gplots)
library(enrichplot)

### read data----
targets.jurkat <- read_tsv("study_design_jurkat_2.txt") 
path <- file.path(targets.jurkat$sample, "abundance.tsv") # set file paths to mapped data

Tx.jurkat <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name")) # annotations and gene symbols 
Tx.jurkat <- as_tibble(Tx.jurkat)
Tx.jurkat <- dplyr::rename(Tx.jurkat, target_id = tx_id)
Tx.jurkat <- dplyr::select(Tx.jurkat, "target_id", "gene_name")
Txi_gene.jurkat <- tximport(path, #imports the data for all samples
                            type = "kallisto", 
                            tx2gene = Tx.jurkat, 
                            txOut = FALSE, # determines whether data represented at transcript or gene level (here read at gene level)
                            countsFromAbundance = "lengthScaledTPM", #transcripts per million
                            ignoreTxVersion = TRUE) 


### preprocessing----

sampleLabels.jurkat <- targets.jurkat$sample
myDGEList.jurkat <- DGEList(Txi_gene.jurkat$counts)
log2.cpm.jurkat <- cpm(myDGEList.jurkat, log=TRUE)  ###log2 normalised counts per million

log2.cpm.jurkat.df <- as_tibble(log2.cpm.jurkat, rownames = "geneID") ###convert as data frame
colnames(log2.cpm.jurkat.df) <- c("geneID", sampleLabels.jurkat)

###tidy data
log2.cpm.jurkat.df.pivot <- pivot_longer(log2.cpm.jurkat.df, 
                                         cols = c(A10, A11, A12, A13, A14, A15),
                                         names_to = "samples", # name of new column
                                         values_to = "expression") # name of new column storing all the data

###plot tidy data of log2 expression - violin plots

p1 <- ggplot(log2.cpm.jurkat.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               colour = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()


###filter data
cpm.jurkat <- cpm(myDGEList.jurkat)
keepers.jurkat <- rowSums(cpm.jurkat>1)>=3 #user defined - depends on studydesign, eliminates lowly expressed genes
myDGEList.jurkat.filtered <- myDGEList.jurkat[keepers.jurkat,]

log2.cpm.jurkat.filtered <- cpm(myDGEList.jurkat.filtered, log=TRUE)
log2.cpm.jurkat.filtered.df <- as_tibble(log2.cpm.jurkat.filtered, rownames = "geneID")
colnames(log2.cpm.jurkat.filtered.df) <- c("geneID", sampleLabels.jurkat)
log2.cpm.jurkat.filtered.df.pivot <- pivot_longer(log2.cpm.jurkat.filtered.df,
                                                  cols = c(A10, A11, A12, A13, A14, A15), 
                                                  names_to = "samples", # name of new column
                                                  values_to = "expression") # name of new column storing all the data

p2 <- ggplot(log2.cpm.jurkat.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

myDGEList.jurkat.filtered.norm <- calcNormFactors(myDGEList.jurkat.filtered, method = "TMM") #TMM normalization
log2.cpm.jurkat.filtered.norm <- cpm(myDGEList.jurkat.filtered.norm, log=TRUE)
log2.cpm.jurkat.filtered.norm.df <- as_tibble(log2.cpm.jurkat.filtered.norm, rownames = "geneID")
colnames(log2.cpm.jurkat.filtered.norm.df) <- c("geneID", sampleLabels.jurkat)
log2.cpm.jurkat.filtered.norm.df.pivot <- pivot_longer(log2.cpm.jurkat.filtered.norm.df,
                                                       cols = c(A10, A11, A12, A13, A14, A15),
                                                       names_to = "samples", # name of new column
                                                       values_to = "expression") # name of new column storing all the data


p3 <- ggplot(log2.cpm.jurkat.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)


### multivariate analysis ----

mydata.jurkat.df <- mutate(log2.cpm.jurkat.filtered.norm.df,
                           jurkat.wt.AVG = (A13 + A11 + A12)/3, 
                           jurkat.a7.AVG = (A10 + A14 + A15)/3,
                           #now make columns comparing each of the averages above 
                           LogFC = (jurkat.a7.AVG - jurkat.wt.AVG)) %>% 
  mutate_if(is.numeric, round, 2)

datatable(mydata.jurkat.df[,c(1,8:10)], 
          extensions = c('KeyTable', "FixedHeader"), 
          filter = 'top',
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         lengthMenu = c("10", "25", "50", "100")))


groupjurkat <- targets.jurkat$phenotype
groupjurkat <- factor(groupjurkat)


### pca jurkat ----

pca.res.jurkat <- prcomp(t(log2.cpm.jurkat.filtered.norm), scale.=F, retx=T)
pc.var.jurkat <- pca.res.jurkat$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per.jurkat <- round(pc.var.jurkat/sum(pc.var.jurkat)*100, 1) 
pca.res.jurkat.df <- as_tibble(pca.res.jurkat$x)
pca.plot.jurkat <- ggplot(pca.res.jurkat.df) +
  aes(x=PC1, y=PC2, label=sampleLabels.jurkat, color = groupjurkat) +
  geom_point(size=4) +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per.jurkat[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per.jurkat[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

pca.plot.jurkat

###diff genes - volcano plot ----

designjurkat <- model.matrix(~0 + groupjurkat)
colnames(designjurkat) <- levels(groupjurkat)

v.DEGList.jurkat.filtered.norm <- voom(myDGEList.jurkat.filtered.norm, designjurkat, plot = TRUE) #models the mean-variance trend
fit.jurkat <- lmFit(v.DEGList.jurkat.filtered.norm, designjurkat) #fits linear model
contrast.matrix.jurkat <- makeContrasts(differences = Jurkat_a7 - Jurkat_wt, 
                                        levels=designjurkat) #contrast matrix

fits.jurkat <- contrasts.fit(fit.jurkat, contrast.matrix.jurkat) #contrasts between the 2 groups
ebFit.jurkat <- eBayes(fits.jurkat)
myTopHits.jurkat <- topTable(ebFit.jurkat, adjust ="BH", coef=1, number=20000, sort.by="logFC")

myTopHits.jurkat.df <- myTopHits.jurkat %>%
  as_tibble(rownames = "geneID")
gt(myTopHits.jurkat.df) #a table, not necessary

datatable(myTopHits.jurkat.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Table 1: DEGs in Jurkat_wt vs a7',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 25, lengthMenu = c("10", "25", "50", "100")))

write.csv(myTopHits.jurkat.df, file = "Jurkat_wt_v_a7_all_genes.csv", quote = FALSE, row.names = FALSE)

signif.genes <- myTopHits.jurkat.df %>%
  filter(adj.P.Val < 0.05 & (logFC <= -1.5 | logFC >= 1.5))

write.csv(signif.genes, file = "signif_DEGs_jurkat_wt_v_a7_1dot5logfc.csv", quote = FALSE, row.names = FALSE)

a7_signif_genes <- signif.genes #variable will  used in venns_diff_expressed_genes_jurkat.R
a7_signif_genes$up_down <- c(NA)
a7_signif_genes$up_down[a7_signif_genes[, "logFC"] < 0] <- "down"
a7_signif_genes$up_down[a7_signif_genes[, "logFC"] >= 0] <- "up"


vplot <- ggplot(myTopHits.jurkat.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  #geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=1) +
  #geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  #geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  annotate("rect", xmin = 1.5, xmax = 13, ymin = -log10(0.05), ymax = 3.2, alpha=.2, fill="#BE684D") +
  annotate("rect", xmin = -1.5, xmax = -13, ymin = -log10(0.05), ymax = 3.2, alpha=.2, fill="#2C467A") +
  annotate(geom = "text", x = -11, y = 1.1, label = "297", size = 7, colour = "#2C467A") + #labels are manually added after running summary(results.jurkat downstream)
  annotate(geom = "text", x = 11, y = 1.1, label = "383", size = 7, colour = "#BE684D") +
  labs(title="Jurkat wt & a7",
       subtitle = "Volcano plot",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw(base_size = 15)
vplot
ggplotly(vplot) #interactive version - can hover and see each gene


results.jurkat <- decideTests(ebFit.jurkat, method="global", adjust.method="BH", p.value=0.05, lfc=1.5)
summary(results.jurkat) #results are added manually to vplot upstream

colnames(v.DEGList.jurkat.filtered.norm$E) <- sampleLabels.jurkat
diffGenes.jurkat <- v.DEGList.jurkat.filtered.norm$E[results.jurkat[,1] !=0,]
diffGenes.jurkat.df <- as_tibble(diffGenes.jurkat, rownames = "geneID")
datatable(diffGenes.jurkat.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Table 1: DEGs in Jurkat wt vs a7',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:7), digits=2)

write.csv(diffGenes.jurkat.df, file = "DEGs_jurkat_wt_v_a7_per_sample_1dot5logfc.csv", quote = FALSE, row.names = FALSE)

# heatmap with up and downregulated genes----

myheatcolours <- rev(brewer.pal(name="PiYG", n=10))
clustRows.jurkat <- hclust(as.dist(1-cor(t(diffGenes.jurkat), method="pearson")), method="complete") #cluster rows by pearson correlation
clustColumns.jurkat <- hclust(as.dist(1-cor(diffGenes.jurkat, method="spearman")), method="complete")
module.assign.jurkat <- cutree(clustRows.jurkat, k=2)
module.colour.jurkat <- rainbow(length(unique(module.assign.jurkat)), start=0.1, end=0.9) 
module.colour.jurkat <- module.colour.jurkat[as.vector(module.assign.jurkat)] 
heatmap.2(diffGenes.jurkat, 
          Rowv=as.dendrogram(clustRows.jurkat), 
          Colv=as.dendrogram(clustColumns.jurkat),
          RowSideColors=module.colour.jurkat,
          col=myheatcolours, scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(8,20),
          main = "jurkat wt vs a7")


# upregulated genes ----
modulePick.jurkat <- 1 
myModule_up.jurkat <- diffGenes.jurkat[names(module.assign.jurkat[module.assign.jurkat %in% modulePick.jurkat]),] 
hrsub_up.jurkat <- hclust(as.dist(1-cor(t(myModule_up.jurkat), method="pearson")), method="complete") 

heatmap.2(myModule_up.jurkat, 
          Rowv=as.dendrogram(hrsub_up.jurkat), 
          Colv = NA, 
          labRow = NA,
          col=myheatcolours, scale="row", 
          density.info="none", trace="none", 
          RowSideColors=module.colour.jurkat[module.assign.jurkat %in% modulePick.jurkat], margins=c(8,20))

# downregulated genes ----
modulePick.jurkat <- 2 
myModule_down.jurkat <- diffGenes.jurkat[names(module.assign.jurkat[module.assign.jurkat %in% modulePick.jurkat]),] 
hrsub_down.jurkat <- hclust(as.dist(1-cor(t(myModule_down.jurkat), method="pearson")), method="complete") 

heatmap.2(myModule_down.jurkat, 
          Rowv=as.dendrogram(hrsub_down.jurkat), 
          Colv = NA, 
          labRow = NA,
          col=myheatcolours, scale="row", 
          density.info="none", trace="none", 
          RowSideColors=module.colour.jurkat[module.assign.jurkat %in% modulePick.jurkat], margins=c(8,20))

# gostplots - go enrichment analysis ----

gost.res_up.jurkat <- gost(rownames(myModule_up.jurkat), organism = "hsapiens", correction_method = "fdr")
gostplot(gost.res_up.jurkat, interactive = T, capped = F) #,title = "GO term analysis in genes upregulated in jurkat_clones")

gost.res_down.jurkat <- gost(rownames(myModule_down.jurkat), organism = "hsapiens", correction_method = "fdr")
gostplot(gost.res_down.jurkat, interactive = T, capped = F) #, title = "GO term analysis in genes upregulated in jurkat_wt")

# GSEA ----
hs_gsea_c2 <- msigdbr(species = "Homo sapiens", # 
                      category = "C2") %>% # choose msigdb collection of interest
  dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols of genes in each signature 


# the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis
mydata.jurkat.df.sub <- dplyr::select(mydata.jurkat.df, geneID, LogFC)
mydata.jurkat.gsea <- mydata.jurkat.df.sub$LogFC
names(mydata.jurkat.gsea) <- as.character(mydata.jurkat.df.sub$geneID)
mydata.jurkat.gsea <- sort(mydata.jurkat.gsea, decreasing = TRUE)

# GSEA function from clusterProfiler
set.seed(1234)
myGSEA.jurkat.res <- GSEA(mydata.jurkat.gsea, TERM2GENE=hs_gsea_c2, verbose=FALSE, seed = TRUE, nPermSimple = 10000, eps = 0)
myGSEA.jurkat.df <- as_tibble(myGSEA.jurkat.res@result)

# interactive table
datatable(myGSEA.jurkat.df[,1:11],  ##add column 11 for gene list
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Signatures enriched in clones',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(3:8), digits=4) #


gseaplot2(myGSEA.jurkat.res, 
          geneSetID = 1, #can choose multiple signatures to overlay in this plot. corresponds to ID in datatable
          pvalue_table = FALSE, 
          title = myGSEA.jurkat.res$Description[1])


myGSEA.jurkat.df <- myGSEA.jurkat.df %>%
  mutate(phenotype = dplyr::case_when(
    NES > 0 ~ "jurkat_a7",
    NES < 0 ~ "jurkat_wt"))

# bubble plot
ggplot(myGSEA.jurkat.df[1:30, ], aes(x=phenotype, y=ID)) + 
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()
