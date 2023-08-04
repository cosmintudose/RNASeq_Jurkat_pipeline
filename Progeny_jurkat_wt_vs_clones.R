#PROGENy analysis Jurkat WT vs clones
#This code follows the tutorial available here: https://bioc.ism.ac.jp/packages/3.14/bioc/vignettes/progeny/inst/doc/progenyBulk.html
#PROGENy is available here: https://saezlab.github.io/progeny/

#This script calculates PROGENy pathway scores for all pathways included in the database
#Input: Normalised expression dataframe from Jurkat_wt_vs_a7_vs_c9.R
#Input: Differential expression analysis output from Jurkat_wt_vs_clones.R
#Input: Study design
#Output: plot for NES heatmap for all samples across all pathways; 
#Output: plot for NES difference for each pathway, comparing WT vs clones 
#Output: pathway responsive genes for each pathway - scatterplots and heatmap only for WNT

library(progeny)
library(pheatmap)
library(tidyverse)
library(ggrepel)
library(gridExtra)
library(gplots)
library(RColorBrewer)

#load data
Normalised_counts <- read_csv("normalised_expression_jurkat_wt_vs_a7_c9.csv")
Experimental_design <- read_tsv("study_design_jurkat.txt")
diff_expr <- read_csv("Jurkat_wt_v_clones_all_genes.csv")


Normalised_counts_matrix <- Normalised_counts %>% 
  dplyr::mutate_if(~ any(is.na(.x)),~ if_else(is.na(.x),0,.x)) %>% 
  tibble::column_to_rownames(var = "geneID") %>% 
  as.matrix()


diff_expr_wt_vs_clones_matrix <- diff_expr %>% 
  dplyr::select(geneID, t) %>% 
  dplyr::filter(!is.na(t)) %>% 
  column_to_rownames(var = "geneID") %>%
  as.matrix()

PathwayActivity_counts <- progeny(Normalised_counts_matrix, scale=TRUE, 
                                  organism="Human", top = 100)
Activity_counts <- as.vector(PathwayActivity_counts)


paletteLength <- 100
myColor <- colorRampPalette(c("blue", "whitesmoke", "red"))(paletteLength)

progenyBreaks <- c(seq(min(Activity_counts), 0, 
                       length.out=ceiling(paletteLength/2) + 1),
                   seq(max(Activity_counts)/paletteLength, 
                       max(Activity_counts), 
                       length.out=floor(paletteLength/2)))

#heatmap for NES for each samples
progeny_hmap <- pheatmap(t(PathwayActivity_counts),fontsize=14, 
                         fontsize_row = 14, fontsize_col = 14, 
                         color=myColor, breaks = progenyBreaks, 
                         main = "PROGENy (100) - Jurkat - WT vs clones", angle_col = 0,
                         treeheight_row = 0,  border_color = NA)

PathwayActivity_zscore <- progeny(diff_expr_wt_vs_clones_matrix, 
                                  scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = TRUE) %>%
  t()
colnames(PathwayActivity_zscore) <- "NES"

PathwayActivity_zscore_df <- as.data.frame(PathwayActivity_zscore) %>% 
  rownames_to_column(var = "Pathway") %>%
  dplyr::arrange(NES) %>%
  dplyr::mutate(Pathway = factor(Pathway))

#NES difference for all pathways comparing WT and clones
pathways_plot <- ggplot(PathwayActivity_zscore_df,aes(x = reorder(Pathway, -NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "blue", high = "red", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal(base_size = 18) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7), limits = c(-4, 4)) +
  theme(axis.title = element_text(face = "bold", size = 18),
        axis.text.x = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 18)) +
  coord_flip() +
  #labs(title = "Jurkat WT vs clones - enriched pathways") +
  xlab("Pathways")

#library(Cairo)
#ggsave(plot = pathways_plot, file = "./progeny_pathways.pdf", 
 #      width = 6.5, height = 5, dpi = 1000, device = cairo_pdf)


prog_matrix <- getModel("Human", top=100) %>% 
  as.data.frame()  %>%
  tibble::rownames_to_column("GeneID")

diff_expr_wt_vs_clones_df <- diff_expr_wt_vs_clones_matrix %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("geneID") %>%
  left_join(diff_expr, by = "geneID") %>%
  select(c("geneID", "logFC")) %>%
  rename(GeneID = geneID)

write.csv(prog_matrix, file = "jurkat_wt_vs_clones_progeny_weights.csv")

scat_plots <- progenyScatter(df = diff_expr_wt_vs_clones_df,
                             weight_matrix = prog_matrix,
                             statName = "logFC", verbose = FALSE)

#saving scatterplots for all pathways
for (pathway in names(scat_plots[[1]])) {
  png(filename = paste0("progeny_", pathway,
                        "_jurkat_wt_vs_clones.png"),
     width = 750, height = 500)
  plot(scat_plots[[1]][[pathway]]) 
  dev.off()
}

#prepare for WNT heatmap
signif_genes <- diff_expr %>%
  filter(adj.P.Val < 0.25)

wnt_genes <- prog_matrix %>%
  filter(WNT != 0 ) %>%
  select(c("GeneID", "WNT")) %>%
  filter(GeneID %in% signif_genes$geneID) %>%
  rename("geneID" = "GeneID")

merged_wnt_signif_genes <- left_join(wnt_genes, signif_genes)

colnames(Normalised_counts) <- c("geneID", "A7_1", "WT_1", "WT_2", "WT_3", "A7_2", "A7_3", "C9_1", "C9_2", "C9_3")

normalised_counts_wnt <- Normalised_counts %>%
  filter(geneID %in% merged_wnt_signif_genes$geneID) %>%
  as.data.frame()
rownames(normalised_counts_wnt) <- normalised_counts_wnt$geneID
normalised_counts_wnt <- normalised_counts_wnt %>%
  select(-geneID) %>%
  as.matrix()

myheatcolours <- rev(brewer.pal(n = 8, name = "RdBu"))
clustRows.jurkat <- hclust(as.dist(1-cor(t(normalised_counts_wnt), method="pearson")), method="complete") #cluster rows by pearson correlation
clustColumns.jurkat <- hclust(as.dist(1-cor(normalised_counts_wnt, method="pearson")), method="complete")
module.assign.jurkat <- wnt_genes %>%
  mutate(module = case_when(WNT < 0 ~ 1, 
                            WNT > 0 ~ 2)) %>%
  select(-WNT) %>%
  deframe()
module.colour.jurkat <- hcl.colors(length(unique(module.assign.jurkat)), palette = "viridis", rev = TRUE) 
module.colour.jurkat <- module.colour.jurkat[as.vector(module.assign.jurkat)] 

df_labels <- as.data.frame(colnames(normalised_counts_wnt))
df_labels$group <- c("A7", "WT", "WT", "WT", "A7", "A7", "C9", "C9", "C9")

labCol <- df_labels$group[match(colnames(normalised_counts_wnt), df_labels$'colnames(normalised_counts_wnt)') ]

#heatmap
pdf("heatmap_WNT_responsive_genes.pdf", width = 7, height = 4)
heatmap.2(normalised_counts_wnt, 
          dendrogram = 'column',
          Rowv=as.dendrogram(clustRows.jurkat), 
          Colv=as.dendrogram(clustColumns.jurkat),
          #RowSideColors=module.colour.jurkat,
          col=myColor, scale='row', srtCol= 360, adjCol = 0.5, #labRow=NA, 
          density.info="none", trace="none",  
          cexRow=1, cexCol=1.4, margins=c(4,10), labCol = labCol,
          main = "WT vs EZH2-KO\nWNT responsive genes", key.title = NA,
          key.xlab = "Expression\nz-score")
dev.off()
  