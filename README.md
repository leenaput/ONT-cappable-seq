# ONT-cappable-seq

We describe here a general overview of the ONT-cappable-seq data analysis pipeline, including the tools and scripts that were used to analyse the data. 

## **Data generation: library preparation and sequencing**

Total RNA was enriched for primary transcripts using an adapted version of the (SMRT)-cappable-seq enrichment protocol. In parallel, control samples were not enriched but subjected to similar environmental conditions. The RNA was reverse transcribed, PCR amplified and barcoded according to Oxford Nanopore Technology cDNA-PCR protocol (SQK-PCS109 combined with SQK-PBK004). The amplified cDNA samples were pooled together in a final library, which was subsequently loaded on a promethION flowcell (R9.4.1) and sequenced with live base-calling and demultiplexing enabled. Fastq files with Q-scores >7 were retained for further analysis. The overall quality of the sequencing run was assessed using NanoComp (v1.11.2). 

## **Data processing**

Note: many of the different steps outlined in this workflow are inspired by the microbepore pipeline that was used for the analysis of prokaryotic Nanopore RNA-seq data (https://felixgrunberger.github.io/microbepore/).

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
