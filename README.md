# **COMP483 Pipeline Project 2026**

## Description
This project employs Snakemake to automate a pipeline that implements tools introduced in COMP483 Computational Biology to analyze sequencing reads. The samples analyzed in this project are from two Human cytomegalovirus patients, 2 and 6 days post-infection, sequenced by [Cheng et al. 2017](https://pubmed.ncbi.nlm.nih.gov/29158406/). First, the pipeline extracts the coding sequence features and quantifies the TPM of the samples using kallisto. The quantified samples are fed into an R script using the sleuth package to determine the differentially expressed genes between the 2 and 6 dpi samples. Next, only reads that map to the HCMV genome are identified using bowtie2, which are then used to create an assembly for each sample using SPAdes. Finally, the longest assembly contigs are aligned using BLAST+ to find the top 5 strains for each sample. The final results of this analysis of the HCMV samples are summarized in a file named PipelineReport.txt.


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

Sample data are provided in this repo in the sample_data folder. These data were generated from the above sample and paired-end fastq files by writing the first 10,000 reads to a new fastq file. Example:

```
head -n 40000 SRR5660030_1.fastq > sample_data/SRR5660030_1.fastq
```


## Running the Pipeline

First, clone this repo to your local system and move into the repo in your terminal:

```
git clone https://github.com/leahbriscoe830/COMP483_Pipeline_Project_2026.git

cd COMP483_Pipeline_Project_2026
```

Before running the pipeline, make sure all sample test data are saved in a folder named *sample_data* in the cloned repo. Once the repo is cloned to your local system and all dependencies are installed, run the pipeline by simply calling Snakemake and the number of cores to use to run the pipeline. For example, to run the pipeline using 2 cores, type either of the following commands into your terminal:

```
snakemake --cores 2

snakemake -c 2
```

A file called PipelineReport.txt will be generated in your local repo. An example report is available called Briscoe_PipelineReport.txt, generated using all reads from the sample data.


## Citations

Bray, N. L., Pimentel, H., Melsted, P., & Pachter, L. (2016). Near-optimal probabilistic RNA-seq quantification. Nature biotechnology, 34(5), 525–527. https://doi.org/10.1038/nbt.3519

Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., and Madden, T.L. (2009). BLAST+: architecture and applications. BMC Bioinformatics, 10, 421.

Cheng, S., Caviness, K., Buehler, J., Smithey, M., Nikolich-Žugich, J., & Goodrum, F. (2017). Transcriptome-wide characterization of human cytomegalovirus in natural infection and experimental latency. Proceedings of the National Academy of Sciences of the United States of America, 114(49), E10586–E10595. https://doi.org/10.1073/pnas.1710522114

Cock, P. J., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., Friedberg, I., Hamelryck, T., Kauff, F., Wilczynski, B., & de Hoon, M. J. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics (Oxford, England), 25(11), 1422–1423. https://doi.org/10.1093/bioinformatics/btp163

Langmead, B., Wilks, C., Antonescu, V., & Charles, R. (2019). Scaling read aligners to hundreds of threads on general-purpose processors. Bioinformatics (Oxford, England), 35(3), 421–432. https://doi.org/10.1093/bioinformatics/bty648

Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J. (2021). Sustainable data analysis with Snakemake. F1000Res 10, 33. https://doi.org/10.12688/f1000research.29032.1

Pimentel, H., Bray, N. L., Puente, S., Melsted, P., & Pachter, L. (2017). Differential analysis of RNA-seq incorporating quantification uncertainty. Nature methods, 14(7), 687–690. https://doi.org/10.1038/nmeth.4324

Prjibelski, A., Antipov, D., Meleshko, D., Lapidus, A., & Korobeynikov, A. (2020). Using SPAdes De Novo Assembler. Current protocols in bioinformatics, 70(1), e102. https://doi.org/10.1002/cpbi.102

Wickham H, François R, Henry L, Müller K, Vaughan D (2026). dplyr: A Grammar of Data Manipulation. R package version 1.2.0, https://dplyr.tidyverse.org
