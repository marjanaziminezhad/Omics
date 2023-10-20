#!/usr/bin bash

#http://bioinfo.icgeb.res.in/gff/gffdownloads/README_GFF_v2.3.txt
#Installation
#-----------------------------------------------------------------------------
#1. Download GFF_v2.3.tar.gz
#wget http://bioinfo.icgeb.res.in/gff/gffdownloads/GFF_v2.3.tar.gz /home/marjan/Documents/Genomanalyse
#2. gunzip GFF_v2.3.tar.gz
#3. tar -vxf GFF_v2.3.tar
#4. ./install.sh	/home/marjan/Documents/Genomanalyse/GFF_v2.3/install		-------  configures the installtion path and install the program
#6. source /home/marjan/Documents/Genomanalyse/GFF_v2.3/install/gff_scripts/gff_profile	[NOTE:Please run this command prior to the execution of GFF-Ex]
#7. Using same terminal move to the directory containing the input files and run GFF-Ex [GFF_INSTALLATION_PATH/gffex -in <gff_file> -db <genome_file>]


#whole human chrms
wget http://ftp.ensembl.org/pub/release-105/gff3/homo_sapiens/Homo_sapiens.GRCh38.105.chr.gff3.gz /home/marjan/Documents/Genomanalyse
wget http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.alt.fa.gz /home/marjan/Documents/Genomanalyse
#whole mouse genom
wget http://ftp.ensembl.org/pub/release-105/gff3/peromyscus_maniculatus_bairdii/Peromyscus_maniculatus_bairdii.HU_Pman_2.1.105.chr.gff3.gz /home/marjan/Documents/Genomanalyse
wget http://ftp.ensembl.org/pub/release-105/fasta/peromyscus_maniculatus_bairdii/dna/Peromyscus_maniculatus_bairdii.HU_Pman_2.1.dna.primary_assembly.1.fa.gz /home/marjan/Documents/Genomanalyse

for i in *.gz #for each file in the zip folder do the following
do
gunzip -d $i
done 
#8.run Gff-Ex. using command "GFF_INSTALLATION_PATH/gffex -in [gff-file] -db [sequencefile]"		GFF_INSTALLATION_PATH/gffex -in example.gtf -db sequence.fasta OR GFF_INSTALLATION_PATH/gffex -in example2.gff -db sequence2.fasta
/home/marjan/Documents/Genomanalyse/GFF_v2.3/install/gff_scripts/gffex -in /home/marjan/Documents/Genomanalyse/Peromyscus_maniculatus_bairdii.HU_Pman_2.1.105.chr.gff3 -db /home/marjan/Documents/Genomanalyse/Peromyscus_maniculatus_bairdii.HU_Pman_2.1.dna.primary_assembly.1.fa
/home/marjan/Documents/Genomanalyse/GFF_v2.3/install/gff_scripts/gffex -in Homo_sapiens.GRCh38.105.chr.gff3 -db Homo_sapiens.GRCh38.dna.alt.fa
#doesnt work 100%