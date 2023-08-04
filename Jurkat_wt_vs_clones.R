###Jurkat wt vs clones
#This script performs differential gene expression analysis and GSEA on Jurkat EZH2-WT vs Jurkat-EZH2-KO
#Input: kallisto abundance for each sample - available on GEO
#Input: study design file available with the code
#Output: volcano plot, GSEA plots
#Output: differentially expressed genes comparing WT with both A7 and C9 EZH2-KO clones
#Output: files for downstream PROGENy analysis

###load packages ----
library(tidyverse)
library(tximport)
library(ensembldb) 
library(EnsDb.Hsapiens.v86)
library(edgeR)
library(matrixStats)
library(limma)
library(clusterProfiler) 
library(msigdbr) 
library(enrichplot)


### read data----
design <- read_tsv("study_design_jurkat.txt") 
path <- file.path("/home/cosmin/rna_seq_clones/raw_data/", design$sample, "abundance.tsv") # set file paths to mapped data

Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name")) %>% #gene symbols 
  as_tibble() %>%
  dplyr::rename(target_id = tx_id) %>%
  dplyr::select("target_id", "gene_name")

Tx_gene <- tximport(path, #imports the data for all samples
                     type = "kallisto",
                     tx2gene = Tx,
                     txOut = FALSE, #data represented at gene level rather than transcript
                     countsFromAbundance = "lengthScaledTPM", #transcripts per million
                     ignoreTxVersion = TRUE) 


### preprocessing----
sample_labels <- design$sample
DGEList <- DGEList(Tx_gene$counts)
log2.cpm <- cpm(DGEList, log=TRUE)  ###log2 normalised counts per million

log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID") ###convert as data frame
colnames(log2.cpm.df) <- c("geneID", sample_labels)

###tidy data
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df,
                                  cols = c(A10, A11, A12, A13, A14, A15, A16, A17, A18),
                                  names_to = "samples", # name of new column
                                  values_to = "expression") # name of new column storing all the data



###filter data
cpm <- cpm(DGEList)
keepers <- rowSums(cpm>1)>=3 #user defined - depends on studydesign, eliminates lowly expressed genes
DGEList.filtered <- DGEList[keepers,]

###normalize data
DGEList.filtered.norm <- calcNormFactors(DGEList.filtered, method = "TMM") #TMM normalization
log2.cpm.filtered.norm <- DGEList.filtered.norm %>% 
  cpm(log=TRUE) %>% 
  as_tibble(rownames = "geneID")

colnames(log2.cpm.filtered.norm) <- c("geneID", sample_labels)


### multivariate analysis ----
data.df <- mutate(log2.cpm.filtered.norm,
                    jurkat.wt.AVG = (A13 + A11 + A12)/3,
                    jurkat.clones.AVG = (A10 + A14 + A15 + A16 + A17 + A18)/6,
                    LogFC = (jurkat.clones.AVG - jurkat.wt.AVG)) %>% 
  mutate_if(is.numeric, round, 2)


group <- design$phenotype
group <- factor(group)


write.csv(data.df[ ,c(1:10)], file = "normalised_expression_jurkat_wt_vs_clones.csv", row.names = FALSE)


### Differential gene expression analysis ----
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

v.DEGList.filtered.norm <- voom(DGEList.filtered.norm, design, plot = TRUE) #models the mean-variance trend
fit <- lmFit(v.DEGList.filtered.norm, design) #fits linear model
contrast.matrix <- makeContrasts(differences = Jurkat_c - Jurkat_wt,
                                 levels=design) #contrast matrix

fits <- contrasts.fit(fit, contrast.matrix) #contrasts between the 2 groups
ebFit <- eBayes(fits)
TopHits <- topTable(ebFit, adjust ="BH", coef=1, number=20000, sort.by="logFC")

TopHits.df <- TopHits %>%
  as_tibble(rownames = "geneID")

####input for PROGENy
write.csv(TopHits.df, file = "Jurkat_wt_v_clones_all_genes.csv", quote = FALSE, row.names = FALSE)


###diff genes - volcano plot ----
results <- decideTests(ebFit, method = "global", adjust.method = "BH", p.value = 0.05, lfc = 1)
summary(results) #results are added manually to vplot upstream

vplot <- ggplot(TopHits.df) +
  aes(y = -log10(adj.P.Val), x = logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  annotate("rect", xmin = 1, xmax = 11, ymin = -log10(0.05), ymax = 2.1, alpha=.2, fill="#BE684D") +
  annotate("rect", xmin = -1, xmax = -11, ymin = -log10(0.05), ymax = 2.1, alpha=.2, fill="#2C467A") +
  annotate(geom = "text", x = -10, y = 1.4, label = "10", size = 7, colour = "#2C467A") + #labels are manually added after running summary(results)
  annotate(geom = "text", x = 10, y = 1.4, label = "12", size = 7, colour = "#BE684D") +
  labs(title="Jurkat wt & clones",
       subtitle = "Volcano plot") +
  theme_bw(base_size = 15)
vplot


colnames(v.DEGList.filtered.norm$E) <- sample_labels
diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")

#saving differentially expressed genes
write.csv(diffGenes.df, file = "DEGs_jurkat_wt_v_clones_per_sample_1logfc.csv", quote = FALSE, row.names = FALSE)

# GSEA ----
hs_gsea_c2 <- msigdbr(species = "Homo sapiens", # 
                      category = "C2") %>% # choose msigdb collection of interest: C2, C6 and H are used in Lefeivre et al
  dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols of genes in each signature 


#the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis
data.df.sub <- dplyr::select(data.df, geneID, LogFC)
data.gsea <- data.df.sub$LogFC
names(data.gsea) <- as.character(data.df.sub$geneID)
data.gsea <- sort(data.gsea, decreasing = TRUE)

#run GSEA
set.seed(1234)
GSEA.res <- GSEA(data.gsea, TERM2GENE=hs_gsea_c2, verbose=FALSE, seed = TRUE, nPermSimple = 10000, eps = 0)
GSEA.df <- as_tibble(GSEA.res@result)

write.csv(GSEA.df, file = "wt_vs_clones_GSEA_C2_jurkat_enriched.csv") #change name depending on which collection is run

#plot for individual signature
gseaplot2(GSEA.res, 
          geneSetID = 1, #change number based on the number of the wanted signature in GSEA.df
          pvalue_table = FALSE, 
          title = GSEA.res$Description[1])

#plots significant GSEA
GSEA.df[1:30, ] %>% #to plot entire df, simly remove the square brackets; currently it's plotting the top 30 enriched sets
  ggplot(aes(x=NES, y= fct_reorder(ID, NES), fill = NES)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(low = "blue", high = "red") +
  labs(y = "ID", title = "C2 - WT vs clones") + #title to be changed depending on which collection is run
  theme_minimal() + 
  theme(axis.title = element_text(size = 24), 
        axis.text.x = element_text(size = 24),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        title = element_text(size = 16), 
        panel.grid.minor.x = element_blank()) 
