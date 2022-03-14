#!/bin/bash

WD=$(pwd)/..
clipped=$WD/mapping_data/clipped
TTS=$WD/boundary_data/TTS
genome=$WD/genome_data

unique_peaks_plus=$(cat $TTS/barcode04/peaks/barcode04.3end.plus.LUZ7.peaks.oracle.narrowPeak.counts.normalized.clustered.csv $TTS/peaks/barcode07/barcode07.3end.plus.LUZ7.peaks.oracle.narrowPeak.counts.normalized.clustered.csv $TTS/barcode10/peaks/barcode10.3end.plus.LUZ7.peaks.oracle.narrowPeak.counts.normalized.clustered.csv | sort -t ',' -k13 -n -u | awk -F ',' 'NR>1 {print $13}')
unique_peaks_minus=$(cat $TTS/barcode04/peaks/barcode04.3end.minus.LUZ7.peaks.oracle.narrowPeak.counts.normalized.clustered.csv $TTS/peaks/barcode07/barcode07.3end.minus.LUZ7.peaks.oracle.narrowPeak.counts.normalized.clustered.csv $TTS/barcode10/peaks/barcode10.3end.minus.LUZ7.peaks.oracle.narrowPeak.counts.normalized.clustered.csv | sort -t ',' -k13 -n -u | awk -F ',' '{NR>1 print $13}')

#tranfer bam alignment files to bed files

for i in $(seq -f %02g 1 12)

do
	mkdir $TTS/barcode$i/TE
	bedtools bamtobed -i $clipped/barcode$i/barcode$i.clipped.sorted.bam | sort -k1,1 -k2,2n | grep "NC_013691.1" > $clipped/barcode$i/barcode$i.LUZ7.clipped.sorted.bed
	grep -w "+" $clipped/barcode$i/barcode$i.LUZ7.clipped.sorted.bed > $clipped/barcode$i/barcode$i.LUZ7.plus.bed
	grep -w "-" $clipped/barcode$i/barcode$i.LUZ7.clipped.sorted.bed > $clipped/barcode$i/barcode$i.LUZ7.minus.bed
done


echo "processing + peaks"

for a in $unique_peaks_plus
do

	#for each peak, extract all the reads in the same orientation that start before it and calculate the coverage for each genomic position

	awk -v TTS="$a" '$2 < TTS-10' $clipped/barcode04/barcode04.LUZ7.plus.bed > $TTS/barcode04/TE/$a.TTS.barcode04.plus.bed
	awk -v TTS="$a" '$2 < TTS-10' $clipped/barcode07/barcode07.LUZ7.plus.bed > $TTS/barcode07/TE/$a.TTS.barcode07.plus.bed
	awk -v TTS="$a" '$2 < TTS-10' $clipped/barcode10/barcode10.LUZ7.plus.bed > $TTS/barcode10/TE/$a.TTS.barcode10.plus.bed

	bedtools genomecov -g $genome/LUZ7.genome -i $TTS/barcode04/TE/$a.TTS.barcode04.plus.bed -d > $TTS/barcode04/TE/$a.TTS.plus.coverage
	bedtools genomecov -g $genome/LUZ7.genome -i $TTS/barcode07/TE/$a.TTS.barcode07.plus.bed -d > $TTS/barcode07/TE/$a.TTS.plus.coverage
	bedtools genomecov -g $genome/LUZ7.genome -i $TTS/barcode10/TE/$a.TTS.barcode10.plus.bed -d > $TTS/barcode10/TE/$a.TTS.plus.coverage

	#calculate the coverage drop across the terminator region for sample t=5min
	echo "processing t5"

	awk -v TTS="$a" '$2 >= TTS-20 && $2 <= TTS+20' $TTS/barcode04/TE/$a.TTS.plus.coverage > $TTS/barcode04/TE/$a.TTS.plus.terminator.region.coverage
	awk -v TTS="$a" '$2 > TTS' $TTS/barcode04/TE/$a.TTS.plus.terminator.region.coverage | awk -v TTS="$a" '{total += $3} END {print TTS, TTS+20, total/NR}' > $TTS/barcode04/TE/$a.TTS.downstream.average.plus.coverage
	awk -v TTS="$a" '$2 < TTS' $TTS/barcode04/TE/$a.TTS.plus.terminator.region.coverage | awk -v TTS="$a" '{total += $3} END {print TTS, TTS-20, total/NR}' > $TTS/barcode04/TE/$a.TTS.upstream.average.plus.coverage
	paste -d " " $TTS/barcode04/TE/$a.TTS.upstream.average.plus.coverage $TTS/barcode04/TE/$a.TTS.downstream.average.plus.coverage | awk '{print $1, $2, $3, $4, $5, $6, $6/$3, 1-$6/$3}' > $TTS/barcode04/TE/$a.plus.coverage.drop


	echo "processing t10"
  	awk -v TTS="$a" '$2 >= TTS-20 && $2 <= TTS+20' $TTS/barcode07/TE/$a.TTS.plus.coverage > $TTS/barcode07/TE/$a.TTS.plus.terminator.region.coverage
        awk -v TTS="$a" '$2 > TTS' $TTS/barcode07/TE/$a.TTS.plus.terminator.region.coverage | awk -v TTS="$a" '{total += $3} END {print TTS, TTS+20, total/NR}' > $TTS/barcode07/TE/$a.TTS.downstream.average.plus.coverage
        awk -v TTS="$a" '$2 < TTS' $TTS/barcode07/TE/$a.TTS.plus.terminator.region.coverage | awk -v TTS="$a" '{total += $3} END {print TTS, TTS-20, total/NR}' > $TTS/barcode07/TE/$a.TTS.upstream.average.plus.coverage
        paste -d " " $TTS/barcode07/TE/$a.TTS.upstream.average.plus.coverage $TTS/barcode07/TE/$a.TTS.downstream.average.plus.coverage | awk '{print $1, $2, $3, $4, $5, $6, $6/$3, 1-$6/$3}' > $TTS/barcode07/TE/$a.plus.coverage.drop	


	echo "processing t20"
	awk -v TTS="$a" '$2 >= TTS-20 && $2 <= TTS+20' $TTS/barcode10/TE/$a.TTS.plus.coverage > $TTS/barcode10/TE/$a.TTS.plus.terminator.region.coverage
        awk -v TTS="$a" '$2 > TTS' $TTS/barcode10/TE/$a.TTS.plus.terminator.region.coverage | awk -v TTS="$a" '{total += $3} END {print TTS, TTS-20, total/NR}' > $TTS/barcode10/TE/$a.TTS.downstream.average.plus.coverage
        awk -v TTS="$a" '$2 < TTS' $TTS/barcode10/TE/$a.TTS.plus.terminator.region.coverage | awk -v TTS="$a" '{total += $3} END {print TTS, TTS+20, total/NR}' > $TTS/barcode10/TE/$a.TTS.upstream.average.plus.coverage
        paste -d " " $TTS/barcode10/TE/$a.TTS.upstream.average.plus.coverage $TTS/barcode10/TE/$a.TTS.downstream.average.plus.coverage | awk '{print $1, $2, $3, $4, $5, $6, $6/$3, 1-$6/$3}' > $TTS/barcode10/TE/$a.plus.coverage.drop

done 

echo "merging TE files"

cat $TTS/barcode04/TE/*.plus.coverage.drop > $TTS/TE_terminators_5minutes_plus.csv
cat $TTS/barcode07/TE/*.plus.coverage.drop > $TTS/TE_terminators_10minutes_plus.csv
cat $TTS/barcode10/TE/*.plus.coverage.drop > $TTS/TE_terminators_20minutes_plus.csv


echo "deleting temp files"

rm $TTS/barcode04/TE/*.plus.coverage.drop $TTS/barcode04/TE/*.plus.coverage $TTS/barcode04/TE/*.plus.bed
rm $TTS/barcode07/TE/*.plus.coverage.drop $TTS/barcode07/TE/*.plus.coverage $TTS/barcode07/TE/*.plus.bed
rm $TTS/barcode10/TE/*.plus.coverage.drop $TTS/barcode10/TE/*.plus.coverage $TTS/barcode10/TE/*.plus.bed

echo "processing - peaks"

for b in $unique_peaks_minus
do

	#for each peak, extract all the reads in the same orientation that start before it and calculate the coverage for each genomic position

	awk -v TTS="$b" '$3 > TSS+10' $clipped/barcode04/barcode04.LUZ7.minus.bed > $TTS/barcode04/TE/$b.TTS.barcode04.minus.bed
	awk -v TTS="$b" '$3 > TTS+10' $clipped/barcode07/barcode07.LUZ7.minus.bed > $TTS/barcode07/TE/$b.TTS.barcode07.minus.bed
	awk -v TTS="$b" '$3 > TTS+10' $clipped/barcode10/barcode10.LUZ7.minus.bed > $TTS/barcode10/TE/$b.TTS.barcode10.minus.bed

	bedtools genomecov -g $genome/LUZ7.genome -i $TTS/barcode04/TE/$b.TTS.barcode04.minus.bed -d > $TTS/barcode04/TE/$b.TTS.minus.coverage
	bedtools genomecov -g $genome/LUZ7.genome -i $TTS/barcode07/TE/$b.TTS.barcode07.minus.bed -d > $TTS/barcode07/TE/$b.TTS.minus.coverage
	bedtools genomecov -g $genome/LUZ7.genome -i $TTS/barcode10/TE/$b.TTS.barcode10.minus.bed -d > $TTS/barcode10/TE/$b.TTS.minus.coverage

	#calculate the coverage drop across the terminator region for sample t=5min
	echo "processing t5"

	awk -v TTS="$b" '$2 >= TTS-20 && $2 <= TTS+20' $TTS/barcode04/TE/$b.TTS.minus.coverage > $TTS/barcode04/TE/$b.TTS.minus.terminator.region.coverage
	awk -v TTS="$b" '$2 < TTS' $TTS/barcode04/TE/$b.TTS.minus.terminator.region.coverage | awk -v TTS="$b" '{total += $3} END {print TTS, TTS-20, total/NR}' > $TTS/barcode04/TE/$b.TTS.downstream.average.minus.coverage
	awk -v TTS="$b" '$2 > TTS' $TTS/barcode04/TE/$b.TTS.minus.terminator.region.coverage | awk -v TTS="$b" '{total += $3} END {print TTS, TTS+20, total/NR}' > $TTS/barcode04/TE/$b.TTS.upstream.average.minus.coverage
	paste -d " " $TTS/barcode04/TE/$b.TTS.upstream.average.minus.coverage $TTS/barcode04/TE/$b.TTS.downstream.average.minus.coverage | awk '{print $1, $2, $3, $4, $5, $6, $6/$3, 1-$6/$3}' > $TTS/barcode04/TE/$b.minus.coverage.drop


	echo "processing t10"
  	awk -v TTS="$b" '$2 >= TTS-20 && $2 <= TTS+20' $TTS/barcode07/TE/$b.TTS.minus.coverage > $TTS/barcode07/TE/$b.TTS.minus.terminator.region.coverage
        awk -v TTS="$b" '$2 < TTS' $TTS/barcode07/TE/$b.TTS.minus.terminator.region.coverage | awk -v TTS="$b" '{total += $3} END {print TTS, TTS-20, total/NR}' > $TTS/barcode07/TE/$b.TTS.downstream.average.minus.coverage
        awk -v TTS="$b" '$2 > TTS' $TTS/barcode07/TE/$b.TTS.minus.terminator.region.coverage | awk -v TTS="$b" '{total += $3} END {print TTS, TTS+20, total/NR}' > $TTS/barcode07/TE/$b.TTS.upstream.average.minus.coverage
        paste -d " " $TTS/barcode07/TE/$b.TTS.upstream.average.minus.coverage $TTS/barcode07/TE/$b.TTS.downstream.average.minus.coverage | awk '{print $1, $2, $3, $4, $5, $6, $6/$3, 1-$6/$3}' > $TTS/barcode07/TE/$b.minus.coverage.drop	


	echo "processing t20"
	awk -v TTS="$b" '$2 >= TTS-20 && $2 <= TTS+20' $TTS/barcode10/TE/$b.TTS.minus.coverage > $TTS/barcode10/TE/$b.TTS.minus.terminator.region.coverage
        awk -v TTS="$b" '$2 > TTS' $TTS/barcode10/TE/$b.TTS.minus.terminator.region.coverage | awk -v TTS="$b" '{total += $3} END {print TTS, TTS-20, total/NR}' > $TTS/barcode10/TE/$b.TTS.downstream.average.minus.coverage
        awk -v TTS="$b" '$2 < TTS' $TTS/barcode10/TE/$b.TTS.minus.terminator.region.coverage | awk -v TTS="$b" '{total += $3} END {print TTS, TTS+20, total/NR}' > $TTS/barcode10/TE/$b.TTS.upstream.average.minus.coverage
        paste -d " " $TTS/barcode10/TE/$b.TTS.upstream.average.minus.coverage $TTS/barcode10/TE/$b.TTS.downstream.average.minus.coverage | awk '{print $1, $2, $3, $4, $5, $6, $6/$3, 1-$6/$3}' > $TTS/barcode10/TE/$b.minus.coverage.drop

done 

echo "merging TE files"

cat $TTS/barcode04/TE/*.minus.coverage.drop > $TTS/TE_terminators_5minutes_minus.csv
cat $TTS/barcode07/TE/*.minus.coverage.drop > $TTS/TE_terminators_10minutes_minus.csv
cat $TTS/barcode10/TE/*.minus.coverage.drop > $TTS/TE_terminators_20minutes_minus.csv


echo "deleting temp files"

rm $TTS/barcode04/TE/*.minus.coverage.drop $TTS/barcode04/TE/*minus.coverage $TTS/barcode04/TE/*.minus.bed
rm $TTS/barcode07/TE/*.minus.coverage.drop $TTS/barcode07/TE/*minus.coverage $TTS/barcode07/TE/*.minus.bed
rm $TTS/barcode10/TE/*.minus.coverage.drop $TTS/barcode10/TE/*minus.coverage $TTS/barcode10/TE/*.minus.bed

