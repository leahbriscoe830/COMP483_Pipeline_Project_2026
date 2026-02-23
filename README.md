# COMP483_Pipeline_Project_2026

## Description
This project uses Snakemake to automate a pipeline to analyze samples from two Human cytomegalovirus patients 2 and 6 days post-infection, sequenced by Cheng et al. 2017. 

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
The sample data were retrieved from NCBI under SRR5660030, SRR5660033, SRR5660044, and SRR5660045. The SRA location address links were saved to a txt file hcmv_samples.txt. The samples were downloaded to the local system using the following commands:

```
wget -i hcmv_samples.txt
fasterq-dump ./SRR5660030
fasterq-dump ./SRR5660033
fasterq-dump ./SRR5660044
fasterq-dump ./SRR5660045
```

## Running the Pipeline
