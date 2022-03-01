# ONT-cappable-seq

We describe here a general overview of the ONT-cappable-seq data analysis pipeline, including the tools and scripts that were used to analyse the data. 

## **Data generation: library preparation and sequencing**

Total RNA was enriched for primary transcripts using an adapted version of the (SMRT)-cappable-seq enrichment protocol. In parallel, control samples were not enriched but subjected to similar environmental conditions. The RNA was reverse transcribed, PCR amplified and barcoded according to Oxford Nanopore Technology cDNA-PCR protocol (SQK-PCS109 combined with SQK-PBK004). The amplified cDNA samples were pooled together in a final library, which was subsequently loaded on a promethION flowcell (R9.4.1) and sequenced with live base-calling and demultiplexing enabled. Fastq files with Q-scores >7 were retained for further analysis. The overall quality of the sequencing run was assessed using NanoComp (v1.11.2). 

## **Data processing**

Note: many of the steps outlined in this workflow are based on the microbepore pipeline that was used for the analysis of prokaryotic Nanopore RNA-seq data (https://felixgrunberger.github.io/microbepore/), with some modifications tailored to our ONT-cappable-seq approach. 

### **I. Data navigation**

Data files were organized in the following way:

ONT-cappable-seq/
 - scripts
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


note on samples: 
- barcode01 = enriched sample of uninfected sample
- barcode02 = control sample of uninfected sample
- barcode04 = enriched sample of 5 minutes post-infection
- barcode05 = control sample of 5 minutes post-infection
- barcode07 = enriched sample of 10 minutes post-infection
- barcode08 = control sample of 10 minutes post-infection
- barcode010 = enriched sample of 20 minutes post-infection
- barcode011 = control sample of 20 minutes post-infection
        
### **II. Raw read processing**

First, multiple fastq files were merged to a single fastq file per sample. The raw reads were processed by Pychopper (v2.5.0) to identify and correctly orient full-length cDNA reads that contain both SSP and VNP primers.

```bash
#!/bin/bash

WD=$(pwd)
input=$WD/fastq_data/raw
pychopped=$WD/fastq_data/pychopped

for i in $(seq -f %02g 1 12)
do

	echo "processing barcode$i";

	cdna_classifier.py \
 	-r $pychopped/barcode$i/report.pdf \
  	-u $pychopped/barcode$i/unclassified.fastq \
  	-w $pychopped/barcode$i/rescued.fastq \
  	$input/barcode$i/barcode$i.fastq \
  	$pychopped/barcode$i/pychopped_barcode$i.fastq 
  
done 
```


Next, the full-length pychopped reads were trimmed with Cutadapt (v2.7) to remove 3' polyA stretches and sequence remnants of the strands-switching primers. 

```bash
#!/bin/bash

WD=$(pwd)
pychopped=$WD/fastq_data/pychopped
cutAAA=$WD/fastq_data/cut_AAA
cutSSP=$WD/fastq_data/cut_SSP

for i in $(seq -f %02g 1 12)
do

	echo "processing barcode$i";
	
	#remove polyA tail from 3' end
	cutadapt -a "A{10}" -j 0 \
	-o $cutAAA/barcode$i/barcode$i.cutadapt.fastq \
	$pychopped/barcode$i/pychopped_barcode$i.fastq 


  	#remove SSP primer from 5' end
	cutadapt -g "TTTCTGTTGGTGCTGATATTGCTGGG" -j 0 \
	-o $cutSSP/barcode$i/barcode$i.cutadapt_SSP.fastq \
	$cutAAA/barcode$i/barcode$i.cutadapt.fastq 

done
```
### **III. Read mapping to reference genomes**

Minimap2 (v2.17) was used to map the trimmed reads to the reference genomes and the RNA spike-in sequence. In our study, we used the reference genome of _Pseudomonas aeruginosa_ phage LUZ7 (NC_013691.1) and bacterial host _P. aeruginosa_ strain US449 (CP091880). The sequences were concatenated into a single fasta file, references.fasta, which is provided on this page.   

```bash
#!/bin/bash

WD=$(pwd)
cutSSP=$WD/fastq_data/cut_SSP
genome=$WD/genome_data
mapped=$WD/mapping_data/mapped
clipped=$WD/mapping_data/clipped

for i in $(seq -f %02g 1 12)
do

	echo "processing barcode$i";
	
	# make output directory to store mapping data
	mkdir $mapped/barcode$i 
	
	# map reads to reference sequences using minimap2
	minimap2 -ax map-ont -k14 --MD \
	$genome/references.fasta \
	$cutSSP/barcode$i/barcode$i.cutadapt_SSP.fastq > $mapped/barcode$i/barcode$i.sam 
	
done
```

Next, samclip (v0.4.0) was used to exclude reads from the alignment that have more than 10 clipped bases on either side:

```bash
#!/bin/bash

WD=$(pwd)
genome=$WD/genome_data
mapped=$WD/mapping_data/mapped
clipped=$WD/mapping_data/clipped

for i in $(seq -f %02g 1 12)
do

	# make output directory to store mapping data
	mkdir $clipped/barcode$i 
	
	#clipping
	samclip --max 10 --ref $genome/references.fasta < $mapped/barcode$i/barcode$i.sam > $clipped/barcode$i/barcode$i.clipped.sam
	
	# convert alignment format and assess read mapping metrics using samtools (v1.9)
  	samtools view -bS $clipped/barcode$i/barcode$i.clipped.sam -o $clipped/barcode$i/barcode$i.clipped.bam
        samtools sort -o $clipped/barcode$i/barcode$i.clipped.sorted.bam $clipped/barcode$i/barcode$i.clipped.bam
	samtools index O $clipped/barcode$i/barcode$i.clipped.sorted.bam
	samtools flagstat $clipped/barcode$i/barcode$i.clipped.sorted.bam > $clipped/barcode$i/flagstat.txt
	samtools idxstats $clipped/barcode$i/barcode$i.clipped.sorted.bam > $clipped/barcode$i/idxstats.txt

done
```

### **IV. Read counting to genomic features**

Next, featureCounts (v2.0.1) was used to assign the reads to the genomic features of the phage and the host. The gtf annotation files used are provided in this repository. 

```bash
#!/bin/bash

WD=$(pwd)
genome=$WD/genome_data
clipped=$WD/mapping_data/clipped
counts=$WD/featurecount_data

for i in $(seq -f %02g 1 12)
do
	mkdir $counts/barcode$i
	
	# assign reads to host features
	featureCounts -L -O -a $genome/US449.gtf \
	$clipped/barcode$i/barcode$i.clipped.sorted.bam \
	-o $counts/barcode$i/US449 \
	-t transcript
	
	# assign reads to phage features
	featureCounts -L -O -a $genome/LUZ7.gtf \
	$clipped/barcode$i/barcode$i.clipped.sorted.bam \
	-o $counts/barcode$i/LUZ7 \
	-t transcript
	
done
```

### **V. Identification of phage transcription start sites**

#### Generation of strand-specific bed files
For transcription start site (TSS) detection, we created strand-specific bed files from the bam files that indicate the number of reads that start at each genomic position of the phage genome. For this we used bedtools (v2.29.2).


```bash
#!/bin/bash

WD=$(pwd)
clipped=$WD/mapping_data/clipped
TSS=$WD/boundary_data/TSS

# get 5end position for all reads in sorted bam file

for i in $(seq -f %02g 1 12)
do

	echo "processing barcode$i"
	mkdir $TSS/barcode$i
	bedtools genomecov -ibam $clipped/barcode$i/barcode$i.clipped.sorted.bam -bga -5 -strand - > $TSS/barcode$i/barcode$i.5end.minus.bedgraph

	bedtools genomecov -ibam $clipped/barcode$i/barcode$i.clipped.sorted.bam -bga -5 -strand + > $TSS/barcode$i/barcode$i.5end.plus.bedgraph

done
```

#### Peak calling
These strand-specific input files were then used to find local maxima for 5' read ends using the termseq-peaks.py script, a previously published peak-calling tool (https://github.com/nichd-bspc/termseq-peaks). This tool generates output files in the narrowPeaks format. After filtering for phage data only, peak coverage information was  added using bedtools intersect. 


```bash
#!/bin/bash

WD=$(pwd)
clipped=$WD/mapping_data/clipped
TSS=$WD/boundary_data/TSS

for i in $(seq -f %02g 1 12)
do

#determine peak positions

	termseq_peaks $TSS/barcode$i/barcode$i.5end.plus.bedgraph $TSS/barcode$i/barcode$i.5end.plus.bedgraph --peaks $TSS/barcode$i.5end.plus.peaks --strand +
	termseq_peaks $TSS/barcode$i/barcode$i.5end.minus.bedgraph $TSS/barcode$i/barcode$i.5end.minus.bedgraph --peaks $TSS/barcode$i.5end.minus.peaks --strand -
	
#add counts
	bedtools intersect -wao \
	-a $TSS/barcode$i/barcode$i.5end.plus.peaks.oracle.narrowPeak \
	-b $TSS/barcode$i/barcode$i.5end.plus.bedgraph \
	> $TSS/barcode$i/barcode$i.5end.plus.peaks.oracle.narrowPeak.counts

	bedtools intersect -wao \
	-a $TSS/barcode$i/barcode$i.5end.minus.peaks.oracle.narrowPeak \
	-b $TSS/barcode$i/barcode$i.5end.minus.bedgraph \
	>  $TSS/barcode$i/barcode$i.5end.minus.peaks.oracle.narrowPeak.counts
	
#because we are mainly focusing on the identification of viral TSS, we extracted the data of LUZ7 and stored it in a seperate file for subsequent analysis 
	grep -w "LUZ7" $TSS/barcode$i/barcode$i.5end.plus.peaks.oracle.narrowPeak.counts > $TSS/barcode$i/barcode$i.5end.plus.LUZ7.peaks.oracle.narrowPeak.counts 
	
	grep -w "LUZ7" $TSS/barcode$i/barcode$i.5end.minus.peaks.oracle.narrowPeak.counts > $TSS/barcode$i/barcode$i.5end.minus.LUZ7.peaks.oracle.narrowPeak.counts 
	

#add RPM values 
	total_mapped=$(samtools view -c -F4 $clipped/barcode$i.clipped.sorted.bam)

	awk '{print $1, $2, $3, $4, $5, $6, $7, $8 ,$9, $10, $11, $12, $13, $14, $15, '$total_mapped', 1000000*$14/'$total_mapped'}' $TSS/barcode$i/barcode$i.5end.plus.LUZ7.peaks.oracle.narrowPeak.counts 
	> $TSS/barcode$i/barcode$i.5end.plus.LUZ7.peaks.oracle.narrowPeak.counts.normalized

	awk '{print $1, $2, $3, $4, $5, $6, $7, $8 ,$9, $10, $11, $12, $13, $14, $15, '$total_mapped', 1000000*$14/'$total_mapped'}'  $TSS/barcode$i/barcode$i.5end.minus.LUZ7.peaks.oracle.narrowPeak.counts 
	> $TSS/barcode$i/barcode$i.5end.minus.LUZ7.peaks.oracle.narrowPeak.counts.normalized


done
```

#### Peak clustering
After generating the peak tables, peaks were clustered in R as described earlier (https://github.com/felixgrunberger/microbepore/blob/master/Rscripts/end5_detection.R), with cov_min = 5 and merge_w = 20. For example:

```R

#load dplyr library

library(dplyr)

#load TSS data
barcode01.plus.5end.peaks.LUZ7.counts <- read.table("ONT-cappable-seq/boundary_data/TSS/barcode01/barcode01.5end.plus.LUZ7.peaks.oracle.narrowPeak.counts.normalized")

barcode01.minus.5end.peaks.LUZ7.counts <- read.table("ONT-cappable-seq/boundary_data/TSS/barcode01/barcode01.5end.minus.LUZ7.peaks.oracle.narrowPeak.counts.normalized")

# analysis of + peaks

	## specify strand
	strand <- "+"

	## cluster peaks within 5bp and remove peaks with less than 5 reads

	barcode01.plus.5end.peaks.LUZ7.counts.clustered <-barcode01.plus.5end.peaks.LUZ7.counts %>% dplyr::rename(chr = V1, start_peak = V2, end_peak = V3, prominence = V5, strand_peak = V6, width = V10, start_cov = V12, end_cov = V13, cov = V14, width_cov = V15, mapped_reads = V16, RPM = V17) %>% dplyr::select(-V4, -V7, -V8, -V9,-V11) %>% group_by(start_peak, end_peak) %>% filter(cov == max(cov)) %>% dplyr::mutate(decision_v = ifelse(strand == "+", min(end_cov), max(end_cov))) %>% filter(end_cov == decision_v) %>% ungroup() %>% arrange(end_cov) %>%  dplyr::mutate(index = lag(end_cov, default = 1) + as.integer(5), index1 = cumsum(ifelse(index >= start_cov, 0, 1))+1) %>% group_by(index1) %>% filter(cov == max(cov), cov >= 5) 

# analysis of - peaks

	## specify strand
	strand <- "-"
	
	## cluster peaks within 5bp and remove peaks with less than 5 reads

	barcode01.minus.5end.peaks.LUZ7.counts.clustered <-barcode01.minus.5end.peaks.LUZ7.counts %>% dplyr::rename(chr = V1, start_peak = V2, end_peak = V3, 	prominence = V5, strand_peak = V6, width = V10, start_cov = V12, end_cov = V13, cov = V14, width_cov = V15, mapped_reads = V16, RPM = V17) %>% 		dplyr::select(-V4, -V7, -V8, -V9,-V11) %>% group_by(start_peak, end_peak) %>% filter(cov == max(cov)) %>% dplyr::mutate(decision_v = ifelse(strand == "+", min(end_cov), max(end_cov))) %>% filter(end_cov == decision_v) %>% ungroup() %>% arrange(end_cov) %>%  dplyr::mutate(index = lag(end_cov, default = 1) + as.integer(5), index1 = cumsum(ifelse(index >= start_cov, 0, 1))+1) %>% group_by(index1) %>% filter(cov == max(cov), cov >= 5) 

# write to .csv file
write.csv(barcode01.plus.5end.peaks.LUZ7.clustered, "ONT-cappable-seq/boundary_data/TSS/barcode01/barcode01.5end.plus.LUZ7.peaks.oracle.narrowPeak.counts.normalized.clustered.csv", row.names = FALSE)

write.csv(barcode01.minus.5end.peaks.LUZ7.clustered, "ONT-cappable-seq/boundary_data/TSS/barcode01/barcode01.5end.minus.LUZ7.peaks.oracle.narrowPeak.counts.normalized.clustered.csv", row.names = FALSE)

```

#### Calculation of enrichment ratio
At the final TSS identification step, we determined the enrichment ratio for each peak position by dividing the RPM value of the peak in the enriched samples by the RPM value of the peak position in the corresponding control sample using custom bash scripts. 

For the peaks on the + strand: 
```bash
bash scripts/enriched_TSS_peaks_plus.sh
```

For the peaks on the - strand: 
```bash
bash scripts/enriched_TSS_peaks_minus.sh
```

The output files 
