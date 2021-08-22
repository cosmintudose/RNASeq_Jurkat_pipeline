# RNASeq_Jurkat_pipeline
RNA-Seq pipeline Jurkat cell-line EZH2-KO experiment: alignment, QC, DGE analysis, GSEA, GO

1. Install bash tools

- Kallisto 
```sudo apt install -y kallisto```
If command doesn't work, check for requirements: https://pachterlab.github.io/kallisto/source

- FastQC
```sudo apt install -y fastqc```

- MultiQC
```sudo apt install -y multiqc```


2. Download human genome assembly (cdna)
```
wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
```
3. Run kallisto_quant.sh 

This quantifies the number of reads in each sample, the abundance.tsv file for each sample will be used in the downstream analysis
The multiqc step (check code) will merge the fastqc files into one .html to inspect the QC of the samples
```
sh kallisto_quant.sh
```
4. Install required R libraries in Rstudio
- note: some of them are available in CRAN, whilst most in Bioconductor
tidyverse
tximport
ensembldb 
EnsDb.Hsapien.v86
edgeR
matrixStats
cowplot
DT
gt
plotly
limma
RColorBrewer
GSEABase 
Biobase 
GSVA 
gprofiler2 
clusterProfiler 
msigdbr 
enrichplotgplots


5. Run Jurkat_wt_vs_a7_vs_c9.R for:
- violin plots comparing normalized and filtered data with unnormalized and unfiltered data
- PCA plot

6. Run Jurkat_wt_vs_a7.R and Jurkat_wt_vs_c9.R *could've been merged into just one file, but I prefer them this way
- contain DGE analysis, GSEA, GO

7. Venns_diif_expressed_genes_jurkat.R - For venns showing overlap between wt_vs_a7 AND wt_vs_c9
