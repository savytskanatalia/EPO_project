Pipeline:
1. Get UMI sequence from separate index file and insert into read name; check read quality, trim polyG, but don't drop/trim by quality itself >> may affect the reads derived from modified RNAs. Example:
	nextflow run UMI_ID.nf -c generic.config --my_files "/mydir/*{index,R1,R2}_001.fastq.gz" -with-docker savytskanatalia/quantefication2 --outdir "results"

Requirement: pull docker container savytskanatalia/quantefication2 to reproduce or install fastp version 0.20.1
nextflow version 23.10.1.5891 (script written for dsl2)

2. We had each sample sequenced in two lanes, so before further processing we merge two lanes per sample. Get list of sample prefixes and merge appropriate pairs. Command used:
List sample prefices
	ls -1 | sed -r 's/_L00[0-9]_R[0-9]_PROC.fastq.gz//g' | uniq > lst_samples.txt
Merge files
	while read i; do echo "${i} started"; zcat trimmed/${i}_L001_R1_PROC.fastq.gz trimmed/${i}_L002_R1_PROC.fastq.gz > polyA_RNA/${i}_merged_R1.fastq; gzip polyA_RNA/${i}_merged_R1.fastq ; zcat trimmed/${i}_L001_R2_PROC.fastq.gz trimmed/${i}_L002_R2_PROC.fastq.gz > polyA_RNA/${i}_merged_R2.fastq; gzip polyA_RNA/${i}_merged_R2.fastq ; done < lst_samples.txt

A - General mapping and quantification of total and ribo-tagged RNAseq with the usual expression analysis
3.a. First round of mapping to the reference M.musculus genome:
Map
nextflow run star_pe_dsl1.nf -c star_pe_dsl1.config --reads "polyA_RNA/*R{1,2}.fastq.gz" --outdir "polyA_RNA/results" --star_index "mm39_105bp/" --threads 5 --read_length 104 --sorted 2 -with-docker savytskanatalia/quantefication
Deduplicate
nextflow run umi_dedup_dsl1.nf -c umi_dedup_dsl1.config --bam_files "/mnt/md0/natalia/EPO/polyA_RNA/results/mapped/*.bam" --paired 2 --outdir "/mnt/md0/natalia/EPO/polyA_RNA/results/" -with-docker savytskanatalia/umitools
3.b. Quantification and DE using 
