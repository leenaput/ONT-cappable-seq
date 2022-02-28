# ONT-cappable-seq

We describe here a general overview of the ONT-cappable-seq data analysis pipeline, including the tools and scripts that were used to analyse the data. 

## **Data generation: library preparation and sequencing**

Total RNA was enriched for primary transcripts using an adapted version of the (SMRT)-cappable-seq enrichment protocol. In parallel, control samples were not enriched but subjected to similar environmental conditions. The RNA was reverse transcribed, PCR amplified and barcoded according to Oxford Nanopore Technology cDNA-PCR protocol (SQK-PCS109 combined with SQK-PBK004). The amplified cDNA samples were pooled together in a final library, which was subsequently loaded on a promethION flowcell (R9.4.1) and sequenced with live base-calling and demultiplexing enabled. Fastq files with Q-scores >7 were retained for further analysis. The overall quality of the sequencing run was assessed using NanoComp (v1.11.2). 

## **Data processing**

Note: many of the steps outlined in this workflow are based on the microbepore pipeline that was used for the analysis of prokaryotic Nanopore RNA-seq data (https://felixgrunberger.github.io/microbepore/), with some modifications tailored to our ONT-cappable-seq approach. 

### **Data navigation**

Data files were organized in the following way:

ONT-cappable-seq/
 - genome_data
 - fastq_data/
 	- raw
	- pychopped
	- cut_AAA
	- cut_SSP
 - mapped_data/
	- mapped
	- clipped
 - featurecount_data
 - boundary_data/
 	- TSS 	 
	- TTS

        
### **Raw read processing**

First, multiple fastq files were merged to a single fastq file per sample. The raw reads were processed by Pychopper (v2.5.0) to identify and correctly orient full-length cDNA reads that contain both SSP and VNP primers.

```bash
#!/bin/bash

for i in $(seq -f %02g 1 12)
do

	echo "processing barcode$i";

	cdna_classifier.py \
 	-r ONT-cappable-seq/fastq_data/pychopped/barcode$i/report.pdf \
  	-u ONT-cappable-seq/fastq_data/pychopped/barcode$i/unclassified.fastq \
  	-w ONT-cappable-seq/fastq_data/pychopped/barcode$i/rescued.fastq \
  	ONT-cappable-seq/fastq_data/raw/barcode$i/barcode$i.fastq \
  	ONT-cappable-seq/fastq_data/pychopped/barcode$i/pychopped_barcode$i.fastq 
  
done 
```


Next, the full-length pychopped reads were trimmed with Cutadapt (v2.7) to remove 3' polyA stretches and sequence remnants of the strands-switching primers. 

```bash
#!/bin/bash

for i in $(seq -f %02g 1 12)
do

	echo "processing barcode$i";
	
	#remove polyA tail from 3' end
	cutadapt -a "A{10}" -j 0 \
	-o ONT-cappable-seq/fastq_data/cut_AAA/barcode$i/barcode$i.cutadapt.fastq \
	ONT-cappable-seq/fastq_data/pychopped/barcode$i/pychopped_barcode$i.fastq 


  	#remove SSP primer from 5' end
	cutadapt -g "TTTCTGTTGGTGCTGATATTGCTGGG" -j 0 \
	-o ONT-cappable-seq/fastq_data/cut_SSP/barcode$i/barcode$i.cutadapt_SSP.fastq \
	ONT-cappable-seq/fastq_data/cut_AAA/barcode$i/barcode$i.cutadapt.fastq 

done
```
### **Read mapping to reference genomes**

Minimap2 (v2.17) was used to map the trimmed reads to the reference genomes and the RNA spike-in sequence. In our study, we used the reference genome of _Pseudomonas aeruginosa_ phage LUZ7 (NC_013691.1) and bacterial host _P. aeruginosa_ strain US449 (CP091880). The sequences were concatenated into a single fasta file, references.fasta, which is provided on this page.   

```bash
#!/bin/bash

for i in $(seq -f %02g 1 12)
do

	echo "processing barcode$i";
	
	# make output directory to store mapping data
	mkdir ONT-cappable-seq/mapped_data/mapped/barcode$i 
	
	# map reads to reference sequences using minimap2
	minimap2 -ax map-ont -k14 --MD \
	ONT-cappable-seq/genome_data/references.fasta \
	ONT-cappable-seq/fastq_data/cutadaptSSP/barcode$i/barcode$i.cutadapt_SSP.fastq > ONT-cappable-seq/mapped_data/mapped/barcode$i/barcode$i.sam 
	
done
```

Next, samclip (v0.4.0) was used to exclude reads from the alignment that have more than 10 clipped bases on either side:

```bash
#!/bin/bash

for i in $(seq -f %02g 1 12)
do

	samclip --max 10 --ref ONT-cappable-seq/genome_data/references.fasta < ONT-cappable-seq/mapped_data/mapped/barcode$i/barcode$i.sam > ONT-cappable-seq/mapped_data/clipped/barcode$i/barcode$i.clipped.sam
	
	# convert alignment format and assess read mapping metrics using samtools (v1.9)
  	samtools view -bS ONT-cappable-seq/mapped_data/clipped/barcode$i/barcode$i.clipped.sam -o ONT-cappable-seq/mapped_data/clipped/barcode$i/barcode$i.clipped.bam
        samtools sort -o ONT-cappable-seq/mapped_data/clipped/barcode$i/barcode$i.clipped.sorted.bam ONT-cappable-seq/mapped_data/clipped/barcode$i/barcode$i.clipped.bam
	samtools index ONT-cappable-seq/mapped_data/clipped/barcode$i/barcode$i.clipped.sorted.bam
	samtools flagstat ONT-cappable-seq/mapped_data/clipped/barcode$i/barcode$i.clipped.sorted.bam > ONT-cappable-seq/mapped_data/clipped/barcode$i/flagstat.txt
	samtools idxstats ONT-cappable-seq/mapped_data/clipped/barcode$i/barcode$i.clipped.sorted.bam > ONT-cappable-seq/mapped_data/clipped/barcode$i/idxstats.txt

done
```

### **Read counting to genomic features**

Next, featureCounts (v2.0.1) was used to assign the reads to the genomic features of the phage and the host. The gtf annotation files used are provided in this repository. 

```bash
#!/bin/bash

for i in $(seq -f %02g 1 12)
do
	mkdir ONT-cappable-seq/featurecount_data/barcode$i
	featureCounts -L -O -a ONT-cappable-seq/genome_data/US449.gtf \
	ONT-cappable-seq/mapping_data/clipped/barcode$i/barcode$i.clipped.sorted.bam \
	-o ONT-cappable-seq/featurecount_data/barcode$i/US449 \
	-t transcript
	
	featureCounts -L -O -a ONT-cappable-seq/genome_data/LUZ7.gtf \
	ONT-cappable-seq/mapping_data/clipped/barcode$i/barcode$i.clipped.sorted.bam \
	-o ONT-cappable-seq/featurecount_data/barcode$i/LUZ7 \
	-t transcript
	
done
```

### **Identification of phage transcription start sites**

For transcription start site (TSS) detection, we created strand-specific bed files from the bam files that indicate the number of reads that start at each genomic position of the phage genome. For this we used bedtools (v2.29.2).


```bash
#!/bin/bash

# get 5end position for all reads in sorted bam file

for i in $(seq -f %02g 1 12)
do

	echo "processing barcode$i"
	mkdir ONT-cappable-seq/boundary_data/TSS/barcode$i
	bedtools genomecov -ibam ONT-cappable-seq/mapped_data/clipped/barcode$i/barcode$i.clipped.sorted.bam -bga -5 -strand - > ONT-cappable-seq/boundary_data/TSS/barcode$i/barcode$i.5end.minus.bedgraph

	bedtools genomecov -ibam ONT-cappable-seq/mapped_data/clipped/barcode$i/barcode$i.clipped.sorted.bam -bga -5 -strand + > ONT-cappable-seq/boundary_data/TSS/barcode$i/barcode$i.5end.plus.bedgraph

done
```

These strand-specific input files were then used to find local maxima for 5' read ends using the termseq-peaks.py script, a previously published peak-calling tool (https://github.com/nichd-bspc/termseq-peaks). This tool generates output files in the narrowPeaks format. After filtering for phage data only, peak coverage information was  added using bedtools intersect. 


```bash
#!/bin/bash

for i in $(seq -f %02g 1 12)
do

#determine peak positions

	termseq_peaks ONT-cappable-seq/boundary_data/TSS/barcode$i/barcode$i.5end.plus.bedgraph ONT-cappable-seq/boundary_data/TSS/barcode$i/barcode$i.5end.plus.bedgraph --peaks ONT-cappable-seq/boundary_data/TSS/barcode$i.5end.plus.peaks --strand +
	termseq_peaks ONT-cappable-seq/boundary_data/TSS/barcode$i/barcode$i.5end.minus.bedgraph ONT-cappable-seq/boundary_data/TSS/barcode$i/barcode$i.5end.minus.bedgraph --peaks ONT-cappable-seq/boundary_data/TSS/barcode$i.5end.minus.peaks --strand -
	
#add counts
	bedtools intersect -wao \
	-a ONT-cappable-seq/boundary_data/TSS/barcode$i/barcode$i.5end.plus.peaks.oracle.narrowPeak \
	-b ONT-cappable-seq/boundary_data/TSS/barcode$i/barcode$i.5end.plus.bedgraph \
	>  ONT-cappable-seq/boundary_data/TSS/barcode$i/barcode$i.5end.plus.peaks.oracle.narrowPeak.counts

	bedtools intersect -wao \
	-a ONT-cappable-seq/boundary_data/TSS/barcode$i/barcode$i.5end.minus.peaks.oracle.narrowPeak \
	-b ONT-cappable-seq/boundary_data/TSS/barcode$i/barcode$i.5end.minus.bedgraph \
	>  ONT-cappable-seq/boundary_data/TSS/barcode$i/barcode$i.5end.minus.peaks.oracle.narrowPeak.counts
	
#because we are mainly focusing on the identification of viral TSS, we extracted the data of LUZ7 and stored it in a seperate file for subsequent analysis 
	grep -w "LUZ7" ONT-cappable-seq/boundary_data/TSS/barcode$i/barcode$i.5end.plus.peaks.oracle.narrowPeak.counts > ONT-cappable-seq/boundary_data/TSS/barcode$i/barcode$i.5end.plus.LUZ7.peaks.oracle.narrowPeak.counts 
	
	grep -w "LUZ7" ONT-cappable-seq/boundary_data/TSS/barcode$i/barcode$i.5end.minus.peaks.oracle.narrowPeak.counts > ONT-cappable-seq/boundary_data/TSS/barcode$i/barcode$i.5end.minus.LUZ7.peaks.oracle.narrowPeak.counts 
	

#add RPM values 
	total_mapped=$(samtools view -c -F4 ONT-cappable-seq/mapped_data/clipped/barcode$i.clipped.sorted.bam)

	awk '{print $1, $2, $3, $4, $5, $6, $7, $8 ,$9, $10, $11, $12, $13, $14, $15, '$total_mapped', 1000000*$14/'$total_mapped'}' ONT-cappable-seq/boundary_data/TSS/barcode$i/barcode$i.5end.plus.LUZ7.peaks.oracle.narrowPeak.counts 
	> ONT-cappable-seq/boundary_data/TSS/barcode$i/barcode$i.5end.plus.LUZ7.peaks.oracle.narrowPeak.counts.normalized

	awk '{print $1, $2, $3, $4, $5, $6, $7, $8 ,$9, $10, $11, $12, $13, $14, $15, '$total_mapped', 1000000*$14/'$total_mapped'}' ONT-cappable-seq/boundary_data/TSS/barcode$i/barcode$i.5end.minus.LUZ7.peaks.oracle.narrowPeak.counts 
	> ONT-cappable-seq/boundary_data/TSS/barcode$i/barcode$i.5end.minus.LUZ7.peaks.oracle.narrowPeak.counts.normalized


done
```

After generating the peak tables, peaks were clustered in R as described earlier (https://github.com/felixgrunberger/microbepore/blob/master/Rscripts/end5_detection.R), with cov_min = 5 and merge_w = 20  For example:

```R

#load dplyr library

library(dplyr)

#load TSS data
barcode01.plus.5end.peaks.LUZ7.counts <- read.table("ONT-cappable-seq/boundary_data/TSS/barcode01/barcode01.5end.plus.LUZ7.peaks.oracle.narrowPeak.counts.normalized")

barcode01.minus.5end.peaks.LUZ7.counts <- read.table("ONT-cappable-seq/boundary_data/TSS/barcode01/barcode01.5end.minus.LUZ7.peaks.oracle.narrowPeak.counts.normalized")

# analysis of + peak
## specify strand
strand <- "+"

## cluster peaks within 5bp and remove peaks with less than 5 reads

barcode01.plus.5end.peaks.LUZ7.counts.clustered <-barcode01.plus.5end.peaks.LUZ7.counts %>% dplyr::rename(chr = V1, start_peak = V2, end_peak = V3, prominence = V5, strand_peak = V6, width = V10, start_cov = V12, end_cov = V13, cov = V14, width_cov = V15, mapped_reads = V16, RPM = V17) %>% dplyr::select(-V4, -V7, -V8, -V9,-V11) %>% group_by(start_peak, end_peak) %>% filter(cov == max(cov)) %>% dplyr::mutate(decision_v = ifelse(strand == "+", min(end_cov), max(end_cov))) %>% filter(end_cov == decision_v) %>% ungroup() %>% arrange(end_cov) %>%  dplyr::mutate(index = lag(end_cov, default = 1) + as.integer(5), index1 = cumsum(ifelse(index >= start_cov, 0, 1))+1) %>% group_by(index1) %>% filter(cov == max(cov), cov >= 5) 

# analysis of - peak
## specify strand
strand <- "-"

## cluster peaks within 5bp and remove peaks with less than 5 reads

barcode01.minus.5end.peaks.LUZ7.counts.clustered <-barcode01.minus.5end.peaks.LUZ7.counts %>% dplyr::rename(chr = V1, start_peak = V2, end_peak = V3, prominence = V5, strand_peak = V6, width = V10, start_cov = V12, end_cov = V13, cov = V14, width_cov = V15, mapped_reads = V16, RPM = V17) %>% dplyr::select(-V4, -V7, -V8, -V9,-V11) %>% group_by(start_peak, end_peak) %>% filter(cov == max(cov)) %>% dplyr::mutate(decision_v = ifelse(strand == "+", min(end_cov), max(end_cov))) %>% filter(end_cov == decision_v) %>% ungroup() %>% arrange(end_cov) %>%  dplyr::mutate(index = lag(end_cov, default = 1) + as.integer(5), index1 = cumsum(ifelse(index >= start_cov, 0, 1))+1) %>% group_by(index1) %>% filter(cov == max(cov), cov >= 5) 

# write to .csv file
write.csv(barcode01.plus.5end.peaks.LUZ7.clustered, "ONT-cappable-seq/boundary_data/TSS/barcode01/barcode01.5end.plus.LUZ7.peaks.oracle.narrowPeak.counts.normalized.clustered.csv", row.names = FALSE)

write.csv(barcode01.minus.5end.peaks.LUZ7.clustered, "ONT-cappable-seq/boundary_data/TSS/barcode01/barcode01.5end.minus.LUZ7.peaks.oracle.narrowPeak.counts.normalized.clustered.csv", row.names = FALSE)

```
