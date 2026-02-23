samples = {
    "SRR5660030": "2dpi",
    "SRR5660033": "6dpi",
    "SRR5660044": "2dpi",
    "SRR5660045": "6dpi"
}

rule all:
    input:
        "PipelineReport.txt"

rule pipeline_report:
    input:
        cds="kallisto/hcmv_cds.fasta",
        fdr05="sleuth/fdr05_results.txt",
        sample_name=expand("data/{sample}", sample=samples.keys()), 
        before_bowtie=expand("data/{sample}_1.fastq", sample=samples.keys()),
        after_bowtie=expand("bowtie/results/{sample}_mapped_1.fq.gz", sample=samples.keys()),
        blast_result=expand("blast_plus/results/{sample}.csv", sample=samples.keys())
    output:
        "PipelineReport.txt"
    script:
        "pipeline_report.py"

rule blast_plus:
    input:
        contigs="blast_plus/contigs/{sample}_longest_contig.fasta",
        db1="blast_plus/betaherpesvirinae_db/betaherpesvirinae.ndb",
        db2="blast_plus/betaherpesvirinae_db/betaherpesvirinae.nhr",
        db3="blast_plus/betaherpesvirinae_db/betaherpesvirinae.nin",
        db4="blast_plus/betaherpesvirinae_db/betaherpesvirinae.not",
        db5="blast_plus/betaherpesvirinae_db/betaherpesvirinae.nsq",
        db6="blast_plus/betaherpesvirinae_db/betaherpesvirinae.ntf",
        db7="blast_plus/betaherpesvirinae_db/betaherpesvirinae.nto"
    output:
        "blast_plus/results/{sample}.csv"
    params:
        db="blast_plus/betaherpesvirinae_db/betaherpesvirinae"
    shell:
        "blastn -query {input.contigs} -db {params.db} -out {output} -outfmt '6 sacc pident length qstart qend sstart send bitscore evalue stitle'"

rule blast_plus_db:
    output:
        "blast_plus/betaherpesvirinae_db/betaherpesvirinae.ndb",
        "blast_plus/betaherpesvirinae_db/betaherpesvirinae.nhr",
        "blast_plus/betaherpesvirinae_db/betaherpesvirinae.nin",
        "blast_plus/betaherpesvirinae_db/betaherpesvirinae.not",
        "blast_plus/betaherpesvirinae_db/betaherpesvirinae.nsq",
        "blast_plus/betaherpesvirinae_db/betaherpesvirinae.ntf",
        "blast_plus/betaherpesvirinae_db/betaherpesvirinae.nto"
    run:
        shell("datasets download virus genome taxon Betaherpesvirinae --include genome --filename betaherpesvirinae_dataset.zip ")
        shell("unzip betaherpesvirinae_dataset.zip -d betaherpesvirinae_dataset ")
        shell("makeblastdb -in betaherpesvirinae_dataset/ncbi_dataset/data/genomic.fna -out blast_plus/betaherpesvirinae_db/betaherpesvirinae -title betaherpesvirinae -dbtype nucl")

rule longest_contig:
    input:
        "spades_assemblies/{sample}_assembly/contigs.fasta"
    output:
        "blast_plus/contigs/{sample}_longest_contig.fasta"
    script:
        "longest_contig.py"

rule spades_assembly:
    input:
        fq1="bowtie/results/{sample}_mapped_1.fq.gz",
        fq2="bowtie/results/{sample}_mapped_2.fq.gz"
    output:
        "spades_assemblies/{sample}_assembly/contigs.fasta"
    params:
        out_dir="spades_assemblies/{sample}_assembly/"
    shell:
        "spades.py -k 127 -t 2 --only-assembler -1 {input.fq1} -2 {input.fq2} -o {params.out_dir}"

rule bowtie_map:
    input:
        forward_reads="data/{sample}_1.fastq",
        reverse_reads="data/{sample}_2.fastq",
        in1="bowtie/index/hcmv_bowtie_index.1.bt2",
        in2="bowtie/index/hcmv_bowtie_index.2.bt2",
        in3="bowtie/index/hcmv_bowtie_index.3.bt2",
        in4="bowtie/index/hcmv_bowtie_index.4.bt2",
        inrev1="bowtie/index/hcmv_bowtie_index.rev.1.bt2",
        inrev2="bowtie/index/hcmv_bowtie_index.rev.2.bt2"
    output:
        sam="bowtie/results/{sample}.sam",
        fq1="bowtie/results/{sample}_mapped_1.fq.gz",
        fq2="bowtie/results/{sample}_mapped_2.fq.gz"
    params:
        index="bowtie/index/hcmv_bowtie_index",
        reads="bowtie/results/{sample}_mapped_%.fq.gz"
    shell:
        "bowtie2 --quiet -x {params.index} -1 {input.forward_reads} -2 {input.reverse_reads} -S {output.sam} --al-conc-gz {params.reads}"

rule bowtie_index:
    input:
        "hcmv_dataset/ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna"
    output:
        out1="bowtie/index/hcmv_bowtie_index.1.bt2",
        out2="bowtie/index/hcmv_bowtie_index.2.bt2",
        out3="bowtie/index/hcmv_bowtie_index.3.bt2",
        out4="bowtie/index/hcmv_bowtie_index.4.bt2",
        outrev1="bowtie/index/hcmv_bowtie_index.rev.1.bt2",
        outrev2="bowtie/index/hcmv_bowtie_index.rev.2.bt2"
    params:
        index="bowtie/index/hcmv_bowtie_index"
    shell:
        "bowtie2-build {input} {params.index}"

rule sleuth_diff_genes: 
    input:
        "sleuth/hcmv_sleuth_tab.txt"
    output:
        "sleuth/fdr05_results.txt"
    script:
        "sleuth_diff_genes.R"

rule sleuth_tab:
    input:
	    expand("kallisto/results/{sample}", sample=samples.keys())
    output:
        "sleuth/hcmv_sleuth_tab.txt"
    run:
        with open(output[0], 'w') as table:
            table.write("sample\tcondition\tpath\n")
            for sample, condition in samples.items():
                table.write(f"{sample}\t{condition}\tkallisto/results/{sample}\n")

rule kallisto_quant:
    input:
        forward_reads="data/{sample}_1.fastq",
        reverse_reads="data/{sample}_2.fastq",
        index="kallisto/hcmv_kallisto_index.idx"
    output:
        directory("kallisto/results/{sample}")
    run:
        shell("kallisto quant -i {input.index} -o {output} -b 30 -t 2 {input.forward_reads} {input.reverse_reads}")

rule kallisto_index:
    input:
        "kallisto/hcmv_cds.fasta"
    output:
        "kallisto/hcmv_kallisto_index.idx"
    shell: 
        "kallisto index -i {output} {input}"

rule format_cds:
    input:
        "hcmv_dataset/ncbi_dataset/data/GCF_000845245.1/cds_from_genomic.fna"
    output:
        "kallisto/hcmv_cds.fasta"
    script:
        "kallisto_format_cds.py"

rule hcmv_data:
    output:
        "hcmv_dataset/ncbi_dataset/data/GCF_000845245.1/cds_from_genomic.fna",
        "hcmv_dataset/ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna"
    run:
        shell("datasets download genome accession GCF_000845245.1 --include cds,genome --filename hcmv_dataset.zip ")
        shell("unzip hcmv_dataset.zip -d hcmv_dataset")