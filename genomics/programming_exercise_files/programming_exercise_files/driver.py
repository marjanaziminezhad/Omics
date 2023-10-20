import os


input_dir = 'genomics/genomes'
output_dir = 'genomics/proteins/'
blast_dir = 'marjan/Documents/genomics/blastdb'
ncbi_dir = 'genomics/ncbi_protein/'
file_list1 = [f for f in os.listdir(input_dir) if f.endswith('.fna')]

#prodigal -in file.fasta -out file.faa
#makeblastdb - in protein_file.fasta - dbtype prot - out protein_db

for file in file_list1:

    name_prefix1 = file.split('.fna')[0]
    
    output_file1 = name_prefix1 + '.faa'
    
    command1 = 'prodigal -i ' + input_dir + file + ' -a ' + output_dir + output_file1
    
    
    #print(command1)
   
    #os.system(command1)
   
file_list2 = [f for f in os.listdir(output_dir) if f.endswith('.faa')]

for file in file_list2:
    name_prefix2 = file.split('.faa')[0]
    output_file2 = name_prefix2 + '_db'

    command2 = 'makeblastdb -in ' + output_dir + file + ' -dbtype prot -out ' + blast_dir + output_file2

    #print (command2)
    #os.system(command2)

#blastp -query query_proteins_from_ncbi.fasta -db protein_db db -max_target_seqs 1 -evalue 1E-5 -outfmt "6 qseqid sseqid pident evalue sstart send length" -out blastp.results.txt
file_list3 = [f for f in os.listdir(ncbi_dir) if f.endswith('.fasta')]

for file in file_list3:

    name_prefix3 = file.split('.fasta')[0]
    output_file3 = name_prefix3 + '.txt'

    command3 = 'blastp -query ' + ncbi_dir + file + ' -db ' + blast_dir + 'db' + ' -max_target_seqs 1 -evalue 1E-5 -outfmt "6 qseqid sseqid pident evalue sstart send length" -out ' + output_file3
    #print (command3)
    os.system(command3)