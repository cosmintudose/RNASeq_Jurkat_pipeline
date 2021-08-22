#checks quality of fastq files
fastqc *.gz -t 4
echo "fastqc finished"

#builds index for reference genome
kallisto index -i Homo_sapiens.GRCh38.103.cdna.all.index Homo_sapiens.GRCh38.103.cdna.all.fa
echo "Index built"

#pseudoalignment and output in separate directories
kallisto quant -i Homo_sapiens.GRCh38.103.cdna.all.index --genomebam --gtf Homo_sapiens.GRCh38.103.gtf -o A1 -t 4 A1_1.fq.gz A1_2.fq.gz &> A1.log
kallisto quant -i Homo_sapiens.GRCh38.103.cdna.all.index --genomebam --gtf Homo_sapiens.GRCh38.103.gtf -o A2 -t 4 A2_1.fq.gz A2_2.fq.gz &> A2.log
kallisto quant -i Homo_sapiens.GRCh38.103.cdna.all.index --genomebam --gtf Homo_sapiens.GRCh38.103.gtf -o A3 -t 4 A3_1.fq.gz A3_2.fq.gz &> A3.log
kallisto quant -i Homo_sapiens.GRCh38.103.cdna.all.index --genomebam --gtf Homo_sapiens.GRCh38.103.gtf -o A4 -t 4 A4_1.fq.gz A4_2.fq.gz &> A4.log
kallisto quant -i Homo_sapiens.GRCh38.103.cdna.all.index --genomebam --gtf Homo_sapiens.GRCh38.103.gtf -o A5 -t 4 A5_1.fq.gz A5_2.fq.gz &> A5.log
kallisto quant -i Homo_sapiens.GRCh38.103.cdna.all.index --genomebam --gtf Homo_sapiens.GRCh38.103.gtf -o A6 -t 4 A6_1.fq.gz A6_2.fq.gz &> A6.log
kallisto quant -i Homo_sapiens.GRCh38.103.cdna.all.index --genomebam --gtf Homo_sapiens.GRCh38.103.gtf -o A7 -t 4 A7_1.fq.gz A7_2.fq.gz &> A7.log
kallisto quant -i Homo_sapiens.GRCh38.103.cdna.all.index --genomebam --gtf Homo_sapiens.GRCh38.103.gtf -o A8 -t 4 A8_1.fq.gz A8_2.fq.gz &> A8.log
kallisto quant -i Homo_sapiens.GRCh38.103.cdna.all.index --genomebam --gtf Homo_sapiens.GRCh38.103.gtf -o A9 -t 4 A9_1.fq.gz A9_2.fq.gz &> A9.log
kallisto quant -i Homo_sapiens.GRCh38.103.cdna.all.index --genomebam --gtf Homo_sapiens.GRCh38.103.gtf -o A10 -t 4 A10_1.fq.gz A10_2.fq.gz &> A10.log
kallisto quant -i Homo_sapiens.GRCh38.103.cdna.all.index --genomebam --gtf Homo_sapiens.GRCh38.103.gtf -o A11 -t 4 A11_1.fq.gz A11_2.fq.gz &> A11.log
kallisto quant -i Homo_sapiens.GRCh38.103.cdna.all.index --genomebam --gtf Homo_sapiens.GRCh38.103.gtf -o A12 -t 4 A12_1.fq.gz A12_2.fq.gz &> A12.log
kallisto quant -i Homo_sapiens.GRCh38.103.cdna.all.index --genomebam --gtf Homo_sapiens.GRCh38.103.gtf -o A13 -t 4 A13_1.fq.gz A13_2.fq.gz &> A13.log
kallisto quant -i Homo_sapiens.GRCh38.103.cdna.all.index --genomebam --gtf Homo_sapiens.GRCh38.103.gtf -o A14 -t 4 A14_1.fq.gz A14_2.fq.gz &> A14.log
kallisto quant -i Homo_sapiens.GRCh38.103.cdna.all.index --genomebam --gtf Homo_sapiens.GRCh38.103.gtf -o A15 -t 4 A15_1.fq.gz A15_2.fq.gz &> A15.log
kallisto quant -i Homo_sapiens.GRCh38.103.cdna.all.index --genomebam --gtf Homo_sapiens.GRCh38.103.gtf -o A16 -t 4 A16_1.fq.gz A16_2.fq.gz &> A16.log
kallisto quant -i Homo_sapiens.GRCh38.103.cdna.all.index --genomebam --gtf Homo_sapiens.GRCh38.103.gtf -o A17 -t 4 A17_1.fq.gz A17_2.fq.gz &> A17.log
kallisto quant -i Homo_sapiens.GRCh38.103.cdna.all.index --genomebam --gtf Homo_sapiens.GRCh38.103.gtf -o A18 -t 4 A18_1.fq.gz A18_2.fq.gz &> A18.log


echo "Alignment finished"

#summarizes qc and mapping results in single html
multiqc -d .

#calculates read depth at exon3 of EZH2 - not necessary
echo -n "A1 depth on chr7:148846470-148846598 " > /home/cosmin/rna_seq_clones/depth.txt
samtools depth -r 7:148846470-148846598 /home/cosmin/rna_seq_clones/raw_data/A1/pseudoalignments.bam | sort -k3n,3 | awk '{arr[NR]=$3} END {if (NR%2==1) print arr [(NR+1)/2]; else print(arr[NR/2] + arr[NR/2+1])/2}' >> /home/cosmin/rna_seq_clones/depth.txt
echo -n "A2 depth on chr7:148846470-148846598 " >> /home/cosmin/rna_seq_clones/depth.txt
samtools depth -r 7:148846470-148846598 /home/cosmin/rna_seq_clones/raw_data/A2/pseudoalignments.bam | sort -k3n,3 | awk '{arr[NR]=$3} END {if (NR%2==1) print arr [(NR+1)/2]; else print(arr[NR/2] + arr[NR/2+1])/2}' >> /home/cosmin/rna_seq_clones/depth.txt
echo -n "A3 depth on chr7:148846470-148846598 " >> /home/cosmin/rna_seq_clones/depth.txt
samtools depth -r 7:148846470-148846598 /home/cosmin/rna_seq_clones/raw_data/A3/pseudoalignments.bam | sort -k3n,3 | awk '{arr[NR]=$3} END {if (NR%2==1) print arr [(NR+1)/2]; else print(arr[NR/2] + arr[NR/2+1])/2}' >> /home/cosmin/rna_seq_clones/depth.txt
echo -n "A4 depth on chr7:148846470-148846598 " >> /home/cosmin/rna_seq_clones/depth.txt
samtools depth -r 7:148846470-148846598 /home/cosmin/rna_seq_clones/raw_data/A4/pseudoalignments.bam | sort -k3n,3 | awk '{arr[NR]=$3} END {if (NR%2==1) print arr [(NR+1)/2]; else print(arr[NR/2] + arr[NR/2+1])/2}' >> /home/cosmin/rna_seq_clones/depth.txt
echo -n "A5 depth on chr7:148846470-148846598 " >> /home/cosmin/rna_seq_clones/depth.txt
samtools depth -r 7:148846470-148846598 /home/cosmin/rna_seq_clones/raw_data/A5/pseudoalignments.bam | sort -k3n,3 | awk '{arr[NR]=$3} END {if (NR%2==1) print arr [(NR+1)/2]; else print(arr[NR/2] + arr[NR/2+1])/2}' >> /home/cosmin/rna_seq_clones/depth.txt
echo -n "A6 depth on chr7:148846470-148846598 " >> /home/cosmin/rna_seq_clones/depth.txt
samtools depth -r 7:148846470-148846598 /home/cosmin/rna_seq_clones/raw_data/A6/pseudoalignments.bam | sort -k3n,3 | awk '{arr[NR]=$3} END {if (NR%2==1) print arr [(NR+1)/2]; else print(arr[NR/2] + arr[NR/2+1])/2}' >> /home/cosmin/rna_seq_clones/depth.txt
echo -n "A7 depth on chr7:148846470-148846598 " >> /home/cosmin/rna_seq_clones/depth.txt
samtools depth -r 7:148846470-148846598 /home/cosmin/rna_seq_clones/raw_data/A7/pseudoalignments.bam | sort -k3n,3 | awk '{arr[NR]=$3} END {if (NR%2==1) print arr [(NR+1)/2]; else print(arr[NR/2] + arr[NR/2+1])/2}' >> /home/cosmin/rna_seq_clones/depth.txt
echo -n "A8 depth on chr7:148846470-148846598 " >> /home/cosmin/rna_seq_clones/depth.txt
samtools depth -r 7:148846470-148846598 /home/cosmin/rna_seq_clones/raw_data/A8/pseudoalignments.bam | sort -k3n,3 | awk '{arr[NR]=$3} END {if (NR%2==1) print arr [(NR+1)/2]; else print(arr[NR/2] + arr[NR/2+1])/2}' >> /home/cosmin/rna_seq_clones/depth.txt
echo -n "A9 depth on chr7:148846470-148846598 " >> /home/cosmin/rna_seq_clones/depth.txt
samtools depth -r 7:148846470-148846598 /home/cosmin/rna_seq_clones/raw_data/A9/pseudoalignments.bam | sort -k3n,3 | awk '{arr[NR]=$3} END {if (NR%2==1) print arr [(NR+1)/2]; else print(arr[NR/2] + arr[NR/2+1])/2}' >> /home/cosmin/rna_seq_clones/depth.txt
echo -n "A10 depth on chr7:148846470-148846598 " >> /home/cosmin/rna_seq_clones/depth.txt
samtools depth -r 7:148846470-148846598 /home/cosmin/rna_seq_clones/raw_data/A10/pseudoalignments.bam | sort -k3n,3 | awk '{arr[NR]=$3} END {if (NR%2==1) print arr [(NR+1)/2]; else print(arr[NR/2] + arr[NR/2+1])/2}' >> /home/cosmin/rna_seq_clones/depth.txt
echo -n "A11 depth on chr7:148846470-148846598 " >> /home/cosmin/rna_seq_clones/depth.txt
samtools depth -r 7:148846470-148846598 /home/cosmin/rna_seq_clones/raw_data/A11/pseudoalignments.bam | sort -k3n,3 | awk '{arr[NR]=$3} END {if (NR%2==1) print arr [(NR+1)/2]; else print(arr[NR/2] + arr[NR/2+1])/2}' >> /home/cosmin/rna_seq_clones/depth.txt
echo -n "A12 depth on chr7:148846470-148846598 " >> /home/cosmin/rna_seq_clones/depth.txt
samtools depth -r 7:148846470-148846598 /home/cosmin/rna_seq_clones/raw_data/A12/pseudoalignments.bam | sort -k3n,3 | awk '{arr[NR]=$3} END {if (NR%2==1) print arr [(NR+1)/2]; else print(arr[NR/2] + arr[NR/2+1])/2}' >> /home/cosmin/rna_seq_clones/depth.txt
echo -n "A13 depth on chr7:148846470-148846598 " >> /home/cosmin/rna_seq_clones/depth.txt
samtools depth -r 7:148846470-148846598 /home/cosmin/rna_seq_clones/raw_data/A13/pseudoalignments.bam | sort -k3n,3 | awk '{arr[NR]=$3} END {if (NR%2==1) print arr [(NR+1)/2]; else print(arr[NR/2] + arr[NR/2+1])/2}' >> /home/cosmin/rna_seq_clones/depth.txt
echo -n "A14 depth on chr7:148846470-148846598 " >> /home/cosmin/rna_seq_clones/depth.txt
samtools depth -r 7:148846470-148846598 /home/cosmin/rna_seq_clones/raw_data/A14/pseudoalignments.bam | sort -k3n,3 | awk '{arr[NR]=$3} END {if (NR%2==1) print arr [(NR+1)/2]; else print(arr[NR/2] + arr[NR/2+1])/2}' >> /home/cosmin/rna_seq_clones/depth.txt
echo -n "A15 depth on chr7:148846470-148846598 " >> /home/cosmin/rna_seq_clones/depth.txt
samtools depth -r 7:148846470-148846598 /home/cosmin/rna_seq_clones/raw_data/A15/pseudoalignments.bam | sort -k3n,3 | awk '{arr[NR]=$3} END {if (NR%2==1) print arr [(NR+1)/2]; else print(arr[NR/2] + arr[NR/2+1])/2}' >> /home/cosmin/rna_seq_clones/depth.txt
echo -n "A16 depth on chr7:148846470-148846598 " >> /home/cosmin/rna_seq_clones/depth.txt
samtools depth -r 7:148846470-148846598 /home/cosmin/rna_seq_clones/raw_data/A16/pseudoalignments.bam | sort -k3n,3 | awk '{arr[NR]=$3} END {if (NR%2==1) print arr [(NR+1)/2]; else print(arr[NR/2] + arr[NR/2+1])/2}' >> /home/cosmin/rna_seq_clones/depth.txt
echo -n "A17 depth on chr7:148846470-148846598 " >> /home/cosmin/rna_seq_clones/depth.txt
samtools depth -r 7:148846470-148846598 /home/cosmin/rna_seq_clones/raw_data/A17/pseudoalignments.bam | sort -k3n,3 | awk '{arr[NR]=$3} END {if (NR%2==1) print arr [(NR+1)/2]; else print(arr[NR/2] + arr[NR/2+1])/2}' >> /home/cosmin/rna_seq_clones/depth.txt
echo -n "A18 depth on chr7:148846470-148846598 " >> /home/cosmin/rna_seq_clones/depth.txt
samtools depth -r 7:148846470-148846598 /home/cosmin/rna_seq_clones/raw_data/A18/pseudoalignments.bam | sort -k3n,3 | awk '{arr[NR]=$3} END {if (NR%2==1) print arr [(NR+1)/2]; else print(arr[NR/2] + arr[NR/2+1])/2}' >> /home/cosmin/rna_seq_clones/depth.txt

echo "Depth calculation finished"
echo -ne '\a'