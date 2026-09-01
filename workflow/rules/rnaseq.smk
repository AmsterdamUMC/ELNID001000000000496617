rnaseq_files_df = pd.read_excel(config["data"]["rnaseq_files_xlsx"], header=0)
rnaseq_sampleIDs = rnaseq_files_df["SampleID"].tolist()

rnaseq_fastq_df = pd.concat([
  pd.Series(pd.concat([rnaseq_files_df["Basename"].astype(str) + "_R1",
            rnaseq_files_df["Basename"].astype(str) + "_R2"]), name="old_path"),
  pd.Series(pd.concat([rnaseq_files_df["SampleID"].astype(str) + "_R1",
            rnaseq_files_df["SampleID"].astype(str) + "_R2"]), name="new_path"),
], axis=1)
rnaseq_readfiles = rnaseq_fastq_df["new_path"].tolist()

comparisons_rnaseq = [c for c in comparisons if 'rnaseq' in comparisons[c].get('assays', ['rnaseq', 'olink'])]

# Rules

rule rnaseq_cp_fastq_to_scratch:
  input:
    readfile=lambda w: (rnaseq_fastq_df[rnaseq_fastq_df.new_path == w.readfile].old_path + ".fastq.gz").tolist(),
  output:
    config["scratch_dir"] + "/resources/rnaseq/fastq/{readfile}.fastq.gz",
  conda:
    "../envs/rsync.yaml"
  message:
    "--- Copying {wildcards.readfile} to scratch ---"
  log:
    config["scratch_dir"] + "resources/rnaseq/fastq/cp_{readfile}_to_scratch.log",
  shell:
    'rsync -a "{input.readfile}" "{output}" &> "{log}"'

rule rnaseq_cp_star_reference_to_scratch:
  input:
    config["reference"]["star_index"],
  output:
    directory(config["scratch_dir"] + "/resources/rnaseq/reference"),
  conda:
    "../envs/rsync.yaml"
  message:
    "--- Copying STAR reference index to scratch ---"
  log:
    config["scratch_dir"] + "/resources/rnaseq/cp_star_reference_to_scratch.log",
  shell:
    'rsync -ar "{input}/" "{output}" &> "{log}"'

rule rnaseq_fastqc:
  input:
    config["scratch_dir"] + "/resources/rnaseq/fastq/{readfile}.fastq.gz",
  output:
    fastqcdir=directory("output/rnaseq/fastqc/{readfile}"),
    fastqchtml="output/rnaseq/fastqc/{readfile}/{readfile}_fastqc.html",
    fastqczip="output/rnaseq/fastqc/{readfile}/{readfile}_fastqc.zip",
  conda:
    "../envs/fastqc.yaml"
  log:
    "output/rnaseq/fastqc/{readfile}_fastqc.log",
  message:
    "--- FastQC {wildcards.readfile} ---"
  threads: 8
  shell:
    'fastqc --outdir="{output.fastqcdir}" --threads {threads} "{input}" &> "{log}"'

rule rnaseq_star_pe:
  input:
    reference=config["scratch_dir"] + "/resources/rnaseq/reference",
    fastq_r1=config["scratch_dir"] + "/resources/rnaseq/fastq/{sampleID}_R1.fastq.gz",
    fastq_r2=config["scratch_dir"] + "/resources/rnaseq/fastq/{sampleID}_R2.fastq.gz",
  output:
    sampleid_dir=directory(config["scratch_dir"] + "/output/rnaseq/bam/{sampleID}"),
    bamfile=config["scratch_dir"] + "/output/rnaseq/bam/{sampleID}/{sampleID}_pe_Aligned.sortedByCoord.out.bam",
    sjfile=config["scratch_dir"] + "/output/rnaseq/bam/{sampleID}/{sampleID}_pe_SJ.out.tab",
    logfinalfile=config["scratch_dir"] + "/output/rnaseq/bam/{sampleID}/{sampleID}_pe_Log.final.out",
  conda:
    "../envs/star.yaml"
  log:
    "output/rnaseq/bam/{sampleID}/{sampleID}_STAR.log",
  message:
    "--- STAR: paired-end alignment {wildcards.sampleID} ---"
  threads: 10
  params:
    scratch_dir=config["scratch_dir"],
  shell:
    """
    STAR --runThreadN {threads} \
         --genomeDir "{input.reference}" \
         --readFilesIn "{input.fastq_r1}" "{input.fastq_r2}" \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --limitBAMsortRAM 50000000000 \
         --outFileNamePrefix "{params.scratch_dir}/output/rnaseq/bam/{wildcards.sampleID}/{wildcards.sampleID}_pe_" \
         |& tee -a "{log}"
    """

rule rnaseq_cp_bam_to_base:
  input:
    config["scratch_dir"] + "/output/rnaseq/bam/{sampleID}",
  output:
    sampleid_dir=directory("output/rnaseq/bam/{sampleID}"),
    bamfile="output/rnaseq/bam/{sampleID}/{sampleID}_pe_Aligned.sortedByCoord.out.bam",
    sjfile="output/rnaseq/bam/{sampleID}/{sampleID}_pe_SJ.out.tab",
    logfinalfile="output/rnaseq/bam/{sampleID}/{sampleID}_pe_Log.final.out",
    okfile="output/rnaseq/bam/{sampleID}_rsync.ok"
  conda:
    "../envs/rsync.yaml"
  message:
    "--- Copying {input} to base dir ---"
  log:
    "output/rnaseq/bam/cp_{sampleID}_to_base.log",
  shell:
    """
    rsync -arv "{input}/" "{output.sampleid_dir}" &> "{log}"
    touch "{output.okfile}"
    """

rule rnaseq_samtools_filter:
  input:
    bam="output/rnaseq/bam/{sampleID}/{sampleID}_pe_Aligned.sortedByCoord.out.bam",
  output:
    bam_filtered="output/rnaseq/bam_filtered/{sampleID}/{sampleID}_filtered.bam",
  conda:
    "../envs/samtools.yaml"
  log:
    "output/rnaseq/bam_filtered/{sampleID}/{sampleID}_samtools.log",
  message:
    "--- SAMtools: filtering {wildcards.sampleID} ---"
  threads: 8
  shell:
    # remove unmapped/multi-mapped (260), MAPQ<10, and mitochondrial reads
    """
    samtools view -@ {threads} -S -h -F 260 -q 10 "{input.bam}" \
      | awk '($1 ~ /^@/) || ($3 != "MT") {{ print $0 }}' \
      | samtools view -@ {threads} -b -o {output.bam_filtered} - |& tee -a "{log}"
    """

rule rnaseq_samtools_index:
  input:
    bam="output/rnaseq/bam_filtered/{sampleID}/{sampleID}_filtered.bam",
  output:
    bai="output/rnaseq/bam_filtered/{sampleID}/{sampleID}_filtered.bam.bai",
  conda:
    "../envs/samtools.yaml"
  log:
    "output/rnaseq/bam_filtered/{sampleID}/{sampleID}_index.log",
  message:
    "--- SAMtools: indexing {wildcards.sampleID} ---"
  threads: 8
  shell:
    'samtools index -b -@ {threads} "{input.bam}" |& tee -a "{log}"'

rule rnaseq_featurecounts:
  input:
    bams=expand("output/rnaseq/bam_filtered/{sampleID}/{sampleID}_filtered.bam", sampleID=rnaseq_sampleIDs),
    annotationfile=config["reference"]["gtf"],
  output:
    "output/rnaseq/counts/counts.txt",
  conda:
    "../envs/featurecounts.yaml"
  log:
    "output/rnaseq/counts/featurecounts.log",
  message:
    "--- Subread: featureCounts ---"
  threads: 8
  shell:
    """
    featureCounts -T {threads} \
                  -a "{input.annotationfile}" \
                  -t exon \
                  -g gene_id \
                  -p \
                  -o "{output}" \
                  {input.bams} &> "{log}"
    """

rule rnaseq_multiqc:
  input:
    fastqc=expand("output/rnaseq/fastqc/{readfile}/{readfile}_fastqc.html", readfile=rnaseq_readfiles),
    featurecounts="output/rnaseq/counts/counts.txt",
  output:
    multiqcdir=directory("output/rnaseq/multiqc"),
  conda:
    "../envs/multiqc.yaml"
  log:
    "output/rnaseq/multiqc/multiqc.log",
  message:
    "--- MultiQC ---"
  threads: 24
  shell:
    'multiqc "output/rnaseq" -t {threads} -o "{output.multiqcdir}" &> "{log}"'

rule rnaseq_dds_preparation:
  input:
    counts_txt="output/rnaseq/counts/counts.txt",
    sample_metadata_xlsx=config["data"]["sample_metadata_xlsx"],
  output:
    dds_rds="output/rnaseq/dds/dds.Rds",
  conda:
    "../envs/r-deseq2.yaml"
  log:
    "output/rnaseq/dds/rnaseq_dds_preparation.log",
  message:
    "--- Preparing DESeqDataSet RNAseq ---"
  threads: 1
  resources:
    mem_mb=8000,
  script:
    "../scripts/rnaseq/prepare_dds.R"

rule rnaseq_dds_subset:
  input:
    dds_rds="output/rnaseq/dds/dds.Rds",
  output:
    dds_subset_rds="output/rnaseq/analyses/{comparison}/dds/dds_{comparison}.Rds",
    vst_subset_rds="output/rnaseq/analyses/{comparison}/dds/vst_{comparison}.Rds",
  params:
    filter=lambda wc: config["comparisons"][wc.comparison]["filter"]
  log:
    "output/rnaseq/analyses/{comparison}/dds/rnaseq_dds_{comparison}_subsetting.log",
  conda:
    "../envs/r-deseq2.yaml"
  message:
    "--- Subsetting {wildcards.comparison} RNAseq ---"
  threads: 1
  resources:
    mem_mb=16000,
  script:
    "../scripts/rnaseq/subset_dds.R"

rule rnaseq_differential_expression:
  input:
    dds_rds="output/rnaseq/analyses/{comparison}/dds/dds_{comparison}.Rds",
  output:
    deseq_rds="output/rnaseq/analyses/{comparison}/degs/deseq_{comparison}.Rds",
    degs_csv="output/rnaseq/analyses/{comparison}/degs/degs_{comparison}.csv",
  params:
    comparison=lambda wc: config["comparisons"][wc.comparison],
    deseq2_params=lambda wc: config["analysis_params"]["deseq2"],
  log:
    "output/rnaseq/analyses/{comparison}/degs/rnaseq_differential_expression_{comparison}.log",
  conda:
    "../envs/r-deseq2.yaml"
  message:
    "--- Differential expression {wildcards.comparison} RNAseq ---"
  threads: 1
  resources:
    mem_mb=16000,
  script:
    "../scripts/rnaseq/differential_expression.R"

rule rnaseq_gene_set_enrichment_analysis:
  input:
    msigdb_rds="output/general/gs_msigdb.Rds",
    degs_csv="output/rnaseq/analyses/{comparison}/degs/degs_{comparison}.csv",
  output:
    fgsea_csv="output/rnaseq/analyses/{comparison}/fgsea/fgsea_{comparison}.csv",
  params:
    fgsea_params=lambda wc: config["analysis_params"]["fgsea"],
  log:
    "output/rnaseq/analyses/{comparison}/fgsea/rnaseq_gene_set_enrichment_analysis_{comparison}.log",
  conda:
    "../envs/r-fgsea.yaml"
  message:
    "--- Gene set enrichment analysis {wildcards.comparison} RNAseq ---"
  threads: 1
  resources:
    mem_mb=16000,
  script:
    "../scripts/rnaseq/fgsea.R"
