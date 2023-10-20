#!/bin/bash

'to proceed the prozess we need a loop to go over all read files of samples and pass them over the .sh files.'


set -x 'This command, also known as set -o xtrace, enables the tracing or debugging mode in the script. When this option is set, the shell prints each command it executes along with its expanded arguments. It is useful for understanding the flow of the script and debugging any issues that may arise during execution. The commands are printed to the console with a + sign before them.'

set -e 'This command, also known as set -o errexit or set -o pipefail, enables the option to exit immediately if any command in the script fails. If a command returns a non-zero exit status, indicating an error, the script will terminate. This option helps ensure the script stops execution when an error occurs and can be useful for preventing further issues or data corruption. It is often used in scripts where the successful execution of each command is critical'

for filename in "BSF_0832_HJ52LBBXY_1_TB_HPC7_STAT3_1_RNA_S70379_4_sorted" "BSF_0832_HJ52LBBXY_1_TB_HPC7_STAT3_2_RNA_S70380_5_sorted" "BSF_0832_HJ52LBBXY_1_TB_HPC7_STAT3_3_RNA_S70381_6_sorted" "BSF_0832_HJ52LBBXY_1_TB_HPC7_Y640F_1_RNA_S70385_7_sorted" "BSF_0832_HJ52LBBXY_1_TB_HPC7_Y640F_2_RNA_S70384_8_sorted" "BSF_0832_HJ52LBBXY_1_TB_HPC7_Y640F_3_RNA_S70386_9_sorted" "BSF_0832_HJ52LBBXY_1_TB_HPC7_dsRed_1_RNA_S70382_1_sorted" "BSF_0832_HJ52LBBXY_1_TB_HPC7_dsRed_2_RNA_S70383_2_sorted" "BSF_0832_HJ52LBBXY_1_TB_HPC7_dsRed_3_RNA_S70378_3_sorted"
do

	# 3. run 002 script
	cd hw/002_rawdata/
	./run_bam2fastq.sh ${filename}.bam ../003_bamToFastq.out_fastq/ /apps/anaconda3/envs/omics2023/bin/bamToFastq SE


	# 4. run 003 scripts
	cd ../003_bamToFastq.out_fastq/
	./run_fastqc.sh ${filename}_1.fastq $HOME/hw/004_fastqc.out_html.zip /apps/anaconda3/envs/omics2023/bin/fastqc

	./run_trimmo.sh /apps/anaconda3/envs/omics2023/share/trimmomatic-0.39-2/trimmomatic.jar /apps/anaconda3/envs/omics2023/share/trimmomatic-0.39-2/adapters "TruSeq2-SE.fa" "16" ${filename}_1.fastq "$( pwd )/../005_trimmomatic.out_trimmedFastq.gz" 20 6 20 25


	# 5. run 005 scripts
	cd ../005_trimmomatic.out_trimmedFastq.gz/
	./run_fastqc.sh ${filename}_1_SE.fastq.gz $HOME/hw/006_fastqc.out_html.zip /apps/anaconda3/envs/omics2023/bin/fastqc

	./run_star.sh /proj/courses/transcriptomics/SS23/localmirror/genomes/gencode/mus_musculus/M32/modified/GRCm39.primary_assembly.genome.fa /proj/courses/transcriptomics/SS23/localmirror/genomes/gencode/mus_musculus/M32/modified/gencode.vM32.primary_assembly.annotation.gtf /proj/courses/transcriptomics/SS23/localmirror/index/gencode/mus_musculus/M32/ $( pwd )/../007_star.out_bam /apps/anaconda3/envs/omics2023/bin/STAR /apps/anaconda3/envs/omics2023/bin/samtools ${filename}_1_SE.fastq.gz


	# 6. run 007 scripts
	cd ../007_star.out_bam/
	./run_qualimap.sh /apps/anaconda3/envs/omics2023/bin/qualimap ${filename}_1_SE_Aligned.sortedByCoord.out.bam /proj/courses/transcriptomics/SS23/localmirror/genomes/gencode/mus_musculus/M32/modified/gencode.vM32.primary_assembly.annotation.gtf "$( pwd )/../008_qualimap.out_html" HTML 32G

	./run_featureCounts.sh "/apps/anaconda3/envs/omics2023/bin/featureCounts" ${filename}_1_SE_Aligned.sortedByCoord.out.bam /proj/courses/transcriptomics/SS23/localmirror/genomes/gencode/mus_musculus/M32/modified/gencode.vM32.primary_assembly.annotation.gtf $( pwd )/../009_featurecounts.out_tsv 

	cd 
done
