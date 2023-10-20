!conda install -c bioconda prodigal

prodigal -i input_file.fasta -a protein_file.fasta #-i:  Specify FASTA/Genbank input file (default reads from stdin).
# -a:  Write protein translations to the selected file.
