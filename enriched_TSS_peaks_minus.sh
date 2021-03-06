#!/bin/bash

WD=$(pwd)/..
TSS=$WD/boundary_data/TSS

#extract peak positions	for each sample

for i in $(seq -f %02g 1 12)
do

	awk -F "," '{print $13}' $TSS/barcode$i/barcode$i.5end.minus.LUZ7.peaks.oracle.narrowPeak.counts.normalized.clustered.csv > $TSS/barcode$i/minus_peaks_barcode$i.csv
	mkdir $TSS/barcode$i/tmp

done

#For each timepoint, extract common peak position on + strand for enriched and control sample 

t0_minus=$(grep -Fx -f $TSS/barcode01/minus_peaks_barcode01.csv $TSS/barcode02/minus_peaks_barcode02.csv)
t5_minus=$(grep -Fx -f $TSS/barcode04/minus_peaks_barcode04.csv $TSS/barcode05/minus_peaks_barcode05.csv)
t10_minus=$(grep -Fx -f $TSS/barcode07/minus_peaks_barcode07.csv $TSS/barcode08/minus_peaks_barcode08.csv)
t20_minus=$(grep -Fx -f $TSS/barcode10/minus_peaks_barcode10.csv $TSS/barcode11/minus_peaks_barcode11.csv)


#Extract information of common peaks for timepoint 5 minutes

echo "analysing timepoint 5 minutes"

for a in $t5_minus
do

	awk -F "," '$13 == "'$a'" {print $1, $2, $3, $13, $5, $9, $12}' $TSS/barcode04/barcode04.5end.minus.LUZ7.peaks.oracle.narrowPeak.counts.normalized.clustered.csv > $TSS/barcode04/tmp/$a.minus.barcode04.csv 
	awk -F "," '$13 == "'$a'" {print $1, $2, $3, $13, $5, $9, $12}' $TSS/barcode05/barcode05.5end.minus.LUZ7.peaks.oracle.narrowPeak.counts.normalized.clustered.csv > $TSS/barcode05/tmp/$a.minus.barcode05.csv 

done

	cat $TSS/barcode04/tmp/*.minus.barcode04.csv > $TSS/barcode04/tmp/common_minus_TSS_barcode04.csv
	cat $TSS/barcode05/tmp/*.minus.barcode05.csv > $TSS/barcode05/tmp/common_minus_TSS_barcode05.csv

	paste -d ' ' $TSS/barcode04/tmp/common_minus_TSS_barcode04.csv $TSS/barcode05/tmp/common_minus_TSS_barcode05.csv | awk -F ' ' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $7/$14}' | sort -g -k 2,2 > $TSS/t5_minus_TSS_ratio.csv


#Extract information of common peaks for timepoint 10 minutes

echo "analysing timepoint 10 minutes"

for b in $t10_minus
do

	awk -F "," '$13 == "'$b'" {print $1, $2, $3, $13, $5, $9, $12}' $TSS/barcode07/barcode07.5end.minus.LUZ7.peaks.oracle.narrowPeak.counts.normalized.clustered.csv > $TSS/barcode07/tmp/$b.minus.barcode07.csv
	awk -F "," '$13 == "'$b'" {print $1, $2, $3, $13, $5, $9, $12}' $TSS/barcode08/barcode08.5end.minus.LUZ7.peaks.oracle.narrowPeak.counts.normalized.clustered.csv > $TSS/barcode08.tmp/$b.minus.barcode08.csv

done

	cat $TSS/barcode07/tmp/*.minus.barcode07.csv > $TSS/barcode07/tmp/common_minus_TSS_barcode07.csv
	cat $TSS/barcode08/tmp/*.minus.barcode08.csv > $TSS/barcode08/tmp/common_minus_TSS_barcode08.csv

	paste -d ' ' $TSS/barcode07/tmp/common_minus_TSS_barcode07.csv $TSS/barcode08/tmp/common_minus_TSS_barcode08.csv | awk -F ' ' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $7/$14}' | sort -g -k 2,2 > $TSS/t10_minus_TSS_ratio.csv

#Extract information of common peaks for timepoint 20 minutes

echo "analysing timepoint 20 minutes"

for c in $t20_minus
do

        awk -F "," '$13 == "'$c'" {print $1, $2, $3, $13, $5, $9, $12}' $TSS/barcode10/barcode10.5end.minus.LUZ7.peaks.oracle.narrowPeak.counts.normalized.clustered.csv > $TSS/barcode10/tmp/$b.minus.barcode10.csv
        awk -F "," '$13 == "'$c'" {print $1, $2, $3, $13, $5, $9, $12}' $TSS/barcode11/barcode11.5end.minus.LUZ7.peaks.oracle.narrowPeak.counts.normalized.clustered.csv > $TSS/barcode11/tmp/$b.minus.barcode11.csv

done

        cat $TSS/barcode10/tmp/*.minus.barcode10.csv > $TSS/barcode10/tmp/common_minus_TSS_barcode10.csv
        cat $TSS/barcode11/tmp/*.minus.barcode11.csv > $TSS/barcode11/tmp/common_minus_TSS_barcode11.csv

        paste -d ' ' $TSS/barcode10/tmp/common_minus_TSS_barcode10.csv $TSS/barcode11/tmp/common_minus_TSS_barcode11.csv | awk -F ' ' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $7/$14}' | sort -g -k 2,2 > $TSS/t20_minus_TSS_ratio.csv

#remove temporary files

echo "removing files"

for f in $(seq -f %02g 1 12)
do

	rm -r $TSS/barcode$f/tmp

done
 
