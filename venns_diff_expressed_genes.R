c5_signif_genes
c9_signif_genes

c5_up <- c5_signif_genes$geneID[c5_signif_genes[ , "up_down"] == "up" ]
c9_up <- c9_signif_genes$geneID[c9_signif_genes[ , "up_down"] == "up" ]

c5_down <- c5_signif_genes$geneID[c5_signif_genes[ , "up_down"] == "down" ]
c9_down <- c9_signif_genes$geneID[c9_signif_genes[ , "up_down"] == "down" ]

genes.down <- list(C5 = c5_down, C9 = c9_down)
genes.up <- list(C5 = c5_up, C9 = c9_up)
library(ggVennDiagram)

ggVennDiagram(genes.down, label_size = 7)
ggVennDiagram(genes.up, label_size = 7) +
  scale_fill_gradient(low = "#ff0a54", high = "#fae0e4")


library(VennDiagram)
venn.diagram(genes.down, lwd = 0, cex = 2, cat.cex = 3, print.mode = c("raw", "percent"),
             alpha = c(0.5, 0.5), fill = c("#a9d6e5", "#014f86"), 
             "genes.down.venn.tiff")
venn.diagram(genes.up, lwd = 0, cex = 2, cat.cex = 3, print.mode = c("raw", "percent"),
             alpha = c(0.5, 0.5), fill = c("#ffccd5", "#c9184a"), 
             "genes.up.venn.tiff")
