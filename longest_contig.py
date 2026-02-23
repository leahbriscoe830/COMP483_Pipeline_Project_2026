# Find the longest contig to use in the BLAST search

from Bio.Seq import Seq
from Bio import SeqIO

# Get input ("spades_assemblies/{sample}_assembly/contigs.fasta") and output ("blast_plus/contigs/{sample}_longest_contig.fasta") from snakemake
input = snakemake.input[0]
output = snakemake.output[0]

# Parse records from contigs.fasta for each sample as a list
records = list(SeqIO.parse(input, 'fasta'))

# Find the longest contig
longest_contig = max(records, key=len)

# Write longest contig as a fasta file to the output
with open(output, 'w') as outfile:
    SeqIO.write(longest_contig, outfile, 'fasta')