###Jurkat wt vs a7 vs c9 
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
library(colorspace)

### read data----
targets.jurkat <- read_tsv("study_design_jurkat.txt")
path <- file.path(targets.jurkat$sample, "abundance.tsv") # set file paths to mapped data

Tx.jurkat <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name")) # annotations and gene symbols 
Tx.jurkat <- as_tibble(Tx.jurkat)
Tx.jurkat <- dplyr::rename(Tx.jurkat, target_id = tx_id)
Tx.jurkat <- dplyr::select(Tx.jurkat, "target_id", "gene_name")
Txi_gene.jurkat <- tximport(path, 
                            type = "kallisto", 
                            tx2gene = Tx.jurkat, 
                            txOut = FALSE, # determines whether data represented at transcript or gene level (here read at gene level)
                            countsFromAbundance = "lengthScaledTPM",
                            ignoreTxVersion = TRUE)


### preprocessing----

sampleLabels.jurkat <- targets.jurkat$sample
myDGEList.jurkat <- DGEList(Txi_gene.jurkat$counts)
log2.cpm.jurkat <- cpm(myDGEList.jurkat, log=TRUE)  ###log2 normalised counts per million

log2.cpm.jurkat.df <- as_tibble(log2.cpm.jurkat, rownames = "geneID") ###convert as data frame
colnames(log2.cpm.jurkat.df) <- c("geneID", sampleLabels.jurkat)

###tidy data
log2.cpm.jurkat.df.pivot <- pivot_longer(log2.cpm.jurkat.df, 
                                         cols = c(A10, A11, A12, A13, A14, A15, A16, A17, A18),
                                         names_to = "samples", # name of new column
                                         values_to = "expression") # name of new column storing all the data

###plot tidy data of log2 expression

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
                                                  cols = c(A10, A11, A12, A13, A14, A15, A16, A17, A18), 
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

myDGEList.jurkat.filtered.norm <- calcNormFactors(myDGEList.jurkat.filtered, method = "TMM")
log2.cpm.jurkat.filtered.norm <- cpm(myDGEList.jurkat.filtered.norm, log=TRUE)
log2.cpm.jurkat.filtered.norm.df <- as_tibble(log2.cpm.jurkat.filtered.norm, rownames = "geneID")
colnames(log2.cpm.jurkat.filtered.norm.df) <- c("geneID", sampleLabels.jurkat)
log2.cpm.jurkat.filtered.norm.df.pivot <- pivot_longer(log2.cpm.jurkat.filtered.norm.df,
                                                       cols = c(A10, A11, A12, A13, A14, A15, A16, A17, A18),
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
                           jurkat.c9.AVG = (A16 + A17 + A18)/3,
                           #now make columns comparing each of the averages above 
                           LogFC = (jurkat.c9.AVG - jurkat.wt.AVG)) %>% 
  mutate_if(is.numeric, round, 2)

datatable(mydata.jurkat.df[,c(1,2:10)], 
          extensions = c('KeyTable', "FixedHeader"), 
          filter = 'top',
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         lengthMenu = c("10", "25", "50", "100")))
write.csv(mydata.jurkat.df, file = "jurkat.log2.cpm.norm.wt.a7.c9.csv", quote = FALSE)

groupjurkat <- targets.jurkat$phenotype
groupjurkat <- factor(groupjurkat)


### pca jurkat ----

pca.res.jurkat <- prcomp(t(log2.cpm.jurkat.filtered.norm), scale.=F, retx=T)
pc.var.jurkat <- pca.res.jurkat$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per.jurkat <- round(pc.var.jurkat/sum(pc.var.jurkat)*100, 1) 
pca.res.jurkat.df <- as_tibble(pca.res.jurkat$x)


dark_okabe <- darken(c("#E69F00", "#56B4E9", "#009E73"), amount = 0.2) 

pca.plot.jurkat <- ggplot(pca.res.jurkat.df) +
  aes(x=PC1, y=PC2, label=sampleLabels.jurkat, fill = groupjurkat, colour = groupjurkat, scale = TRUE) +
  geom_point(size=7, stroke = 1, shape = 21) +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per.jurkat[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per.jurkat[3],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw(base_size = 14) +
  scale_fill_manual(labels = c("A7", "C9", "WT"), values = c("#E69F00", "#56B4E9", "#009E73")) +
  scale_colour_manual(labels = c("A7", "C9", "WT"), values = dark_okabe) +
  labs(fill = "Group", colour = "Group")

pca.plot.jurkat

# interactive version
ggplotly(pca.plot.jurkat) 
# 3D plot containing PC3 as well
plot_ly(data = pca.res.jurkat.df, x = pca.res.jurkat.df$PC1, y = pca.res.jurkat.df$PC2, z = pca.res.jurkat.df$PC3, type = "scatter3d", label=sampleLabels.jurkat, color = groupjurkat)
