import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import re
import os

# Function to extract gene IDs from a GFF file
def extract_gene_ids_from_gff(file_path):
    gene_ids = set()
    with open(file_path, 'r') as file:
        for line in file:
            # Assuming gene IDs in GFF files are in a specific format like "ID=Gene12345;"
            match = re.search(r'ID=([^;]+);', line)
            if match:
                gene_id = match.group(1)
                gene_ids.add(gene_id)
    return gene_ids

# Directory paths for the two genomes
elsdenii_directory = "~/cam_proj/06_prokka_e"
hexanoica_directory = "~/cam_proj/06_prokka_h"

# Initialize sets to store gene IDs for each genome
gene_ids_elsdenii = set()
gene_ids_hexanoica = set()

# Iterate over files in the "M. elsdenii" directory
for root, dirs, files in os.walk(elsdenii_directory):
    for file in files:
        if file.endswith(".gff"):
            gff_file_path = os.path.join(root, file)
            gene_ids_elsdenii.update(extract_gene_ids_from_gff(gff_file_path))

# Iterate over files in the "M. hexanoica" directory
for root, dirs, files in os.walk(hexanoica_directory):
    for file in files:
        if file.endswith(".gff"):
            gff_file_path = os.path.join(root, file)
            gene_ids_hexanoica.update(extract_gene_ids_from_gff(gff_file_path))

# Create a Venn diagram
venn2([gene_ids_elsdenii, gene_ids_hexanoica], ('M. elsdenii', 'M. hexanoica'))

# Set labels and title
plt.title("Gene Sets Comparison")
plt.xlabel("Genomes")

# Show the Venn diagram
plt.show()
