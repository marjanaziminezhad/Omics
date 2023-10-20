import os
from muscle_fasta import fetch

input_dir = 'genomics/genomes/'
predicted_proteins = 'genomics/predicted_proteins/'
ncbi_dir = 'genomics/ncbi_protein/'
muscle_dir = 'genomics/muscle_fasta/'
trees_dir = 'genomics/trees/'

file_list1 = [f for f in os.listdir(input_dir) if f.endswith('.fna')]

#prodigal -i input_file.fasta -a protein_file.fasta

for file in file_list1:

    name_prefix1 = file.split('.fna')[0]

    output_file1 = name_prefix1 + '.fasta'

    command1 = 'prodigal -i ' + input_dir + file + ' -a ' + predicted_proteins + output_file1

    #print(command1)
    #os.system(command1)

#makeblastdb - in protein_file.fasta - dbtype prot - out protein_db
file_list2 = [f for f in os.listdir(predicted_proteins) if f.endswith('.fasta')]

for file in file_list2:
    name_prefix2 = file.split('.fasta')[0]
    output_file2 = name_prefix2 + '_db'

    command2 = 'makeblastdb -in ' + predicted_proteins + file + ' -dbtype prot -out ' + output_file2

    #print (command2)
    #os.system(command2)

#blastdbcmd -list blastdb
#blastp -query query_proteins_from_ncbi.fasta -db protein_db db -max_target_seqs 1 -evalue 1E-5 -outfmt "6 qseqid sseqid pident evalue sstart send length" -out blastp.results.txt
file_list3 = [f for f in os.listdir(ncbi_dir) if f.endswith('.fasta')]                                                                   # Query Seq-id, subject seq id,Percentage of identical matches
#i = 1 this is to check if my code loops over all the db files i have
for file in file_list3:
    #print(i)
    #i =i+1
    db_list = [f for f in os.listdir(ncbi_dir) if f.endswith('.phr')] #to have the dbfiles as a list for loop
    #j = 1
    for db in db_list:
        #i = i + 1
        db_name = db.split('.phr')[0] 
        #name_prefix3 = file.split('.fasta')[0]
        #output_file3 = name_prefix3 + '.txt'
        #command3 = 'blastp -query ' + ncbi_dir + file + ' -db ' + ncbi_dir + '/' + db_name + 
        # ' -max_target_seqs 10 -evalue 1E-5 -outfmt "6 qseqid sseqid pident evalue sstart send length" -out ' + 'output_file' + str(i)+ '.txt'
        #print (command3)
        #os.system(command3)

# I need a better overview of my result so I cocat all the results in a file
#os.system('cat *txt > all_text')
# filtering them in POI order
#os.system('grep "CAL35216" all_text > peptidase_Campylobacter_jejuni_subsp.txt')
#os.system('grep "WP_012108828" all_text > cation_proton_antiporter_Campylobacter_hominis.txt')
#os.system('grep "WP_005872928" all_text > signal_peptidase_II_Campylobacter_gracilis.txt')
#os.system('grep "WP_039619574" all_text > DNA-directed_RNA_polymerase_subunit_beta_Campylobacter.txt')
#os.system('grep "WP_044780241" all_text > DNA_topoisomerase_ATP-hydrolyzing_subunit_A_Campylobacter.txt')
#os.system('grep "WP_044794075" all_text > DNA_ligase_Campylobacter.txt')
#print(cd genomics/genomes | cat * fna > fasta_muscle.fasta)
#os.system('cd genomics/genomes |cat *fna > muscle.fasta')
fetch.fetch() # calls the fetch.py to go through text files and build the filtered fasta files
file_list4 = [f for f in os.listdir(muscle_dir) if f.endswith('.txt')]

for file in file_list4:

    name_prefix4 = file.split('.txt')[0]

    output_file4 = name_prefix4 + '.fasta'

    output_muscle = 'aligned_' + name_prefix4
    command4='muscle -in ' + output_file4 + ' -out ' + trees_dir + output_muscle + '.afa'
    print(command4)
    #print(output_file4)
    #os.system(command4)

#need to install the libraries foe iqtree and make a link to in in build
#os.system('wget https://github.com/iqtree/iqtree2/releases/download/v2.2.0/iqtree-2.2.0-Linux.tar.gz')
#os.system('gunzip -d iqtree-2.2.0-Linux.tar.gz')
#os.system('cd /Documents')
#os.system('cmake /home/marjan/Documents/genomics/eigen-3.4.0/')
#os.system('make -j')
#os.system('sudo apt-get install iqtree')

#iqtree -s input.afa
file_list5 = [f for f in os.listdir(trees_dir) if f.endswith('.afa')]

for file in file_list5:
    
 
    name_prefix5 = file.split('.afa')[0]
    output_file5 = name_prefix5
    command5 = 'iqtree -s ' + trees_dir + '/' + file 

    #print (command5)
    #os.system(command5)
