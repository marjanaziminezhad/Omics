#!/usr/bin/env bash

mkdir genomanalyse_bash
cd genomanalyse_bash
#downloading data 
wget https://ftp.ensembl.org/pub/release-108/gff3/homo_sapiens/Homo_sapiens.GRCh38.108.chr.gff3.gz

wget https://ftp.ensembl.org/pub/release-108/gff3/mus_musculus/Mus_musculus.GRCm39.108.gff3.gz 

# annotation in gff3 file $3 -> "exon" & gene - exon -> "intron"
#making new directories for unzip and result of each file
for i in *.gz 
do
mkdir ${i}_unzip
mkdir ${i}_results
mv "$i" "${i}_unzip"
cd "${i}_unzip"
gunzip -d $i #unzip'

grep -v '^#' *.gff3 | sort -t$'\t' -k1,1 -k4,4n -k5,5n > sorted.gff3 #grep -v '^#' < basename $i .gz > nonmatching lines and sorting -t, --field-separator=SEP
              #use SEP instead of non-blank to blank transition & -k, --key=KEYDEF
              #sort via a key; KEYDEF gives location and type
awk '{ if($3 == "chromosome") print $1"\t"$4-1"\t"$5-1 }' sorted.gff3 > TotalChrom.bed #\t new line
awk '{ sum += ($3-$2)} END {print sum}' TotalChrom.bed>chrom_sum_value.txt # -t$'\t' specifies delimiter - TAB for .gff3


awk '{ if($3 == "exon") print $0 }' sorted.gff3 | awk '{print $1"\t"$4"\t"$5} ' | uniq -u > sorted_exon.bed
awk '{ sum += ($3-$2)} END {print sum}' sorted_exon.bed>exon_sum_value.txt

awk '{ if($3 == "gene") print $0 }' sorted.gff3 | awk '{print $1"\t"$4"\t"$5} ' | uniq -u > sorted_gene.bed
awk '{ sum += ($3-$2)} END {print sum}' sorted_gene.bed>gene_sum_value.txt

cat gene_sum_value.txt chrom_sum_value.txt > temp # cat smallervalue biggervlaue > awk 'p{print $0-p}{p=$0}' concatenatedfile
awk 'p{print $0-p}{p=$0}' temp > intergenic_sum_value.txt
rm temp
cat exon_sum_value.txt gene_sum_value.txt >temp1
awk 'p{print $0-p}{p=$0}' temp1 > intron_sum_value.txt
rm temp1
mv *txt ../ 
cd .. 
mv *txt ${i}_results 
cd ${i}_results

echo -e 'Exon_Lenght\tIntron Length\tIntergenic_Length\tTotal_Cromomosome_Length' > final_${i}.txt #header for my output table
paste exon_sum_value.txt intron_sum_value.txt intergenic_sum_value.txt chrom_sum_value.txt >> final_${i}.txt # append all results
awk '{if (NR!=1) {print $1/$4}}' final_${i}.txt > temp1.txt 
awk '{if (NR!=1) {print $2/$4}}' final_${i}.txt > temp2.txt
awk '{if (NR!=1) {print $3/$4}}' final_${i}.txt > temp3.txt 
paste temp1.txt temp2.txt temp3.txt >> final_${i}.txt #rewrite all in finall
rm temp* #remove
cd ..

done
#plot via python matplotlib
chmod +x /home/marjan/Documents/Genomanalyse/bash/plot.py

python /home/marjan/Documents/Genomanalyse/bash/plot.py
