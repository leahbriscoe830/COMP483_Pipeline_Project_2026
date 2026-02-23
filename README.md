# COMP483_Pipeline_Project_2026

## Description
This project uses Snakemake to automate a pipeline implementing tools introduced in COMP483 Computational Biology to analyze sequencing reads. The samples analyzed in this project are from two Human cytomegalovirus patients, 2 and 6 days post-infection, sequenced by Cheng et al. 2017. First, the pipeline extracts the coding sequence features and quantifies the TPM of the samples using kallisto. The quantified samples are fed into an R script using the sleuth package to find the differentially expressed genes between the 2 and 6 dpi samples. Next, only reads that map to the HCMV genome are identified using bowtie2, which are then used to create an assembly for each sample using SPAdes. Finally, the sample assemblies are aligned using blast+ to find the top 5 strains. The final results of this analysis of the HCMV samples are summarized in a file named PipelineReport.txt.

## Dependencies

### Bioinformatics Tools
- NCBI datasets
- kallisto
 -sleuth
- bowtie2
- spades
- blast+

### Python Modules
- os
- csv
- gzip
- Bio.Seq
- Bio.SeqIO

### R Libraries
- sleuth
- dplyr

## Data
The sample data were retrieved from NCBI under the accession numbers SRR5660030, SRR5660033, SRR5660044, and SRR5660045. The SRA location address links were saved to a txt file hcmv_samples.txt. The samples were downloaded to the local system using the following commands:

```
wget -i hcmv_samples.txt
fasterq-dump ./SRR5660030
fasterq-dump ./SRR5660033
fasterq-dump ./SRR5660044
fasterq-dump ./SRR5660045
```
Sample data are provided in this repo in the sample_data folder. These data were generated from the above paired-end fastq files by writing the first 10,000 reads to a new fastq file:

```
head -n 40000 data/SRR5660030_1.fastq > sample_data/SRR5660030_1.fastq
```

## Running the Pipeline
