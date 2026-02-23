# **COMP483 Pipeline Project 2026**

## Description
This project uses Snakemake to automate a pipeline implementing tools introduced in COMP483 Computational Biology to analyze sequencing reads. The samples analyzed in this project are from two Human cytomegalovirus patients, 2 and 6 days post-infection, sequenced by Cheng et al. 2017. First, the pipeline extracts the coding sequence features and quantifies the TPM of the samples using kallisto. The quantified samples are fed into an R script using the sleuth package to find the differentially expressed genes between the 2 and 6 dpi samples. Next, only reads that map to the HCMV genome are identified using bowtie2, which are then used to create an assembly for each sample using SPAdes. Finally, the assembly contigs are aligned using BLAST+ to find the top 5 strains for each sample. The final results of this analysis of the HCMV samples are summarized in a file named PipelineReport.txt.

## Dependencies

### Bioinformatics Tools
- [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- [NCBI datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/)
- [kallisto](https://pachterlab.github.io/kallisto/download)
- [bowtie2](https://github.com/BenLangmead/bowtie2)
- [SPAdes](https://github.com/ablab/spades)
- [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata)

### Python Modules
- [Biopython (Bio.Seq and Bio.SeqIO)](https://biopython.org/)

### R Libraries
- [sleuth](https://pachterlab.github.io/sleuth/download)
- [dplyr](https://dplyr.tidyverse.org/index.html)

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

First, clone this repo to your local system and move into the repo in your terminal

```
git clone https://github.com/leahbriscoe830/COMP483_Pipeline_Project_2026.git
cd COMP483_Pipeline_Project_2026
```

Before running the pipeline, make sure all sample test data are saved in a folder named sample_data in the cloned repo. Once the repo is cloned to your local system and all dependencies are installed, run the pipeline by simply calling Snakemake and the number of cores to use to run the pipeline. For example, to run the pipeline using 2 cores, type either of the following commands into your terminal:

```
snakemake --cores 2
snakemake -c 2
```

A file called PipelineReport.txt will be generated in your local repo. An example report is available called Briscoe_PipelineReport.txt, generated using all reads from the sample data.

## Citations
