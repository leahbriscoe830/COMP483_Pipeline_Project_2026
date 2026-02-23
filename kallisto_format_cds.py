# Change the names of the sequence headers to just the RefSeq protein_id

from Bio.Seq import Seq
from Bio import SeqIO

# Get input ("hcmv_dataset/ncbi_dataset/data/GCF_000845245.1/cds_from_genomic.fna") and output ("kallisto/hcmv_cds.fasta") from snakemake
input = snakemake.input[0]
output = snakemake.output[0]

with open(input, 'r') as infile, open(output, 'w') as outfile:

    # Iterate through the cds fasta
    for line in infile:

        # Only modify headers
        if line.startswith(">"):

            # Extract protein_id
            parts = line.split("[protein_id=")

            # Split the string at the ] to get the protein id
            protein_id = parts[1].split("]")[0]

            # Write the new header to the modified fasta file
            outfile.write(">"+str(protein_id)+"\n")

        # Write the sequence to the modified fasta as normal
        else:
            outfile.write(line) 