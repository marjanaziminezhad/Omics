# FHWN Tulln - Practice - Workflow Bio Data Science

## Workflow Overview

The project consists of several distinct workflows, each serving a specific purpose. Here's an overview of the workflows and their respective steps:

### Raw Data Preparation

1. **Align Reads to Reference Genome**: Align raw sequencing reads to a reference genome using a tool such as STAR.

2. **Index BAM Files**: Create an index of the resulting BAM files to facilitate data retrieval.

3. **Check Alignment**: Ensure the alignment quality and correctness using various quality control tools.

4. **bamToFastq - FASTQ**: Convert BAM files back to FASTQ format if needed.

### Data Quality Assessment

5. **STAR Aligner**: Utilize the STAR aligner to align reads to a reference genome.

6. **FastQC**: Assess the quality of sequencing reads with FastQC.

7. **Trimmomatic**: Perform read trimming and quality filtering using Trimmomatic.

8. **FastQC (Post-Trimming)**: Re-run FastQC to evaluate the quality of reads after trimming.

9. **Samtools Index**: Index the processed BAM files for further analysis.

### Quantification and Differential Expression

10. **Raw Reads - BAM**: Prepare and organize the raw sequencing reads for further analysis.

11. **Qualimap**: Assess the quality of read mapping using Qualimap.

12. **Quantification**: Quantify gene expression using tools like featureCounts.

13. **Differential Expression (DESeq2)**: Perform differential expression analysis using DESeq2.

### Functional Analysis

14. **Functional Analysis (Cluster Profiler)**: Conduct functional analysis using Cluster Profiler to gain insights into biological functions and pathways.

