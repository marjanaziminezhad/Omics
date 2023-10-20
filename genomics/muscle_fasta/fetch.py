import os


def fetch():


    muscle_dir = 'genomics/muscle_fasta/'


    file_list = [f for f in os.listdir(muscle_dir) if f.endswith('.txt')]
    # print(file_list4)


    for file in file_list:
        name_prefix = file.split('.txt')[0]

        output_file = name_prefix + '.fasta'
        # Open the input text file

        with open(os.path.join(muscle_dir, file), "r") as text_file:
            # Loop through each line of the text file
            for line in text_file:
                # Get the second word of the line
                search_term = line.split()[1]
                # print(search_term)
                # Open the FASTA file from predicted_proteind dir.
                f = 'genomics/muscle_fasta/all_predicted.fasta'
                with open(f, "r") as fasta_predicted_protein:
                    # Loop through each line of the FASTA file
                    # for fasta_line in fasta_predicted_protein:
                    # Check if the line contains the search term= 2nd word of text file
                    command4 = 'sed -n ' + '\'/' + search_term + \
                        ' /,/\*/p\' ' + f + ' >> ' + output_file
                    os.system(command4)
                    # print(command4)
