# Use the produced results from the pipeline to generate a pipline report

import csv
import os
import gzip
from Bio.Seq import Seq
from Bio import SeqIO

# Gather all inputs needed to write to the report
cds = snakemake.input.cds
fdr05 = snakemake.input.fdr05
sample_path = list(snakemake.input.sample_name)
before_bowtie = list(snakemake.input.before_bowtie)
after_bowtie = list(snakemake.input.after_bowtie)
blast_result = list(snakemake.input.blast_result)
output = snakemake.output[0]

with open(output, 'w') as outfile:

    ### Write the count of the CDS to the pipeline report file ###

    # Parse CDS fasta
    cds_records = list(SeqIO.parse(cds, 'fasta'))

    # Write the number of cds (records) to the report
    outfile.write("The HCMV genome (GCF_000845245.1) has "+ str(len(cds_records)) +" CDS.\n\n")


    ### Write the details for each significant transcript from Sleuth (FDR < 0.05) ###

    # Open the txt output file and read in the table
    with open(fdr05, 'r') as sleuth_results:
        sig_transcripts = sleuth_results.read()

    # Write the table to the report
    outfile.write(sig_transcripts)


    ### Write the number of reads in each sample before and after the Bowtie2 mapping ###

    # Use zip() to iterate through each sample reads
    for name, before, after in zip(sample_path, before_bowtie, after_bowtie):

        # Get the name of the sample
        sample_name =  os.path.basename(str(name))
        
        # Parse the paired end samples before mapping with bowtie
        sample_records = list(SeqIO.parse(before, 'fastq'))

        # Unzip the fq.gz files and parse the samples after mapping with bowtie
        with gzip.open(after, 'rt') as handle:
            mapped_records = list(SeqIO.parse(handle, 'fastq'))

        # Write the before and after reads for each sample
        outfile.write("\nSample "+ str(sample_name) +" had "+ str(len(sample_records)) +" read pairs before and "+ str(len(mapped_records)) +" read pairs after Bowtie2 filtering.")

    # Formatting
    outfile.write("\n")


    ### Write the two header rows, followed by the top 5 hits, tab-delimit each item ###

    # Use zip() to iterate through each sample blast
    for name, blast in zip(sample_path, blast_result):

        # Get the name of the sample
        sample_name =  os.path.basename(str(name))

        # Read in the blast results csv file
        with open(blast, mode='r', newline='', encoding='utf-8') as blast_csv:
            
            # Write the header indicating the statitstic
            outfile.write("\n"+ str(sample_name) +":\n sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n")
            
            # Set up csv reader, delimit by tab
            csv_reader = csv.reader(blast_csv, delimiter='\t')

            # Iterate through each of the top 5 hits
            row_count = 0
            for row in csv_reader:
                if row_count >= 5:
                    break

                # Read in the row from csv reader and delimt by tab, convert to string
                line = "\t".join(row) 
                
                # Write the row to the report
                outfile.write(line + '\n')
                
                # Increment row count
                row_count += 1