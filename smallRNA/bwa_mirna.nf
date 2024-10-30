// bwa mapping to mature and hairpin RNAs
// part of the miRNA/small RNA analysis
// need to provide bwa index
// also save separately unmapped reads as .fastq
// samtools view -b -f 4 alignments.bam > tmp.bam
// samtools fastq {} > {.}.fastq
// map with bwa, pipe output into samtools to make bam, sort by coordinate; then subset to unmapped only which can go into separate process for .fastq of unmapped generation;
// can be separate processes in a workflow for hairpin and mature
// example for bwa+nf dsl2: https://stackoverflow.com/questions/72967648/bwa-fail-to-load-index-using-nextflow
params.fastq_short = "/some/path/*.fastq.gz"
params.fastq_long = "/some/path/*.fastq.gz"
params.bwa_idx_mature = "data/genome.fa.{,amb,ann,bwt,pac,sa}"
params.bwa_idx_hairpin = "data/genome.fa.{,amb,ann,bwt,pac,sa}"
params.outdir = "results"


process BWA_MAP_SHORT {
    publishDir "${params.outdir}/mapped/below31bpReads", mode: 'move', pattern: "*_sorted.bam"

    input:
    tuple val(sample_id), path(fastq)
    path bwa_idx_mature

    output:
    path("*_sorted.bam")
    tuple val(sample_id), path("${sample_id}_unmapped.bam"), emit: bam_un_short


    script:
    def idxbase = bwa_idx_mature[0].baseName
    """
    bwa mem -t 3 "${idxbase}" "${fastq}" | samtools view -@ 3 -b - | samtools sort -@ 3 - > "${sample_id}_sorted.bam"
    samtools view  -@ 3 -b -f 4  "${sample_id}_sorted.bam" >  "${sample_id}_unmapped.bam"

    """
}
process BWA_MAP_LONG {
    publishDir "${params.outdir}/mapped/longerReads", mode: 'move', pattern: "*_sorted.bam"

    input:
    tuple val(sample_id), path(fastq)
    path bwa_idx_mature

    output:
    path("*_sorted.bam")
    tuple val(sample_id), path("${sample_id}_unmapped.bam"), emit: bam_un_long


    script:
    def idxbase = bwa_idx_mature[0].baseName
    """
    bwa mem -t 3 "${idxbase}" "${fastq}" | samtools view -@ 3 -b - | samtools sort -@ 3 - > "${sample_id}_sorted.bam"
    samtools view  -@ 3 -b -f 4  "${sample_id}_sorted.bam" >  "${sample_id}_unmapped.bam"

    """
}

process SAMTOOLS_FQ_s {

    publishDir "${params.outdir}/unmapped/below31bpReads", mode: 'move', pattern: "*.fastq.gz"

    input:
    tuple val(sample_id), path(bam)

    output:
    path("*.fastq.gz")

    script:

    """
    samtools fastq -@ 3 "${bam}" > "${sample_id}.fastq.gz"


    """
}
process SAMTOOLS_FQ_l {
    publishDir "${params.outdir}/unmapped/longerReads", mode: 'move', pattern: "*.fastq.gz"

    input:
    tuple val(sample_id), path(bam)

    output:
    path("*fastq.gz")

    script:

    """
    samtools fastq -@ 3 "${bam}" > "${sample_id}.fastq.gz"


    """
}

workflow {

    fastq_data_s = channel.fromPath(params.fastq_short).map { file -> tuple(file.baseName, file)}
    bwa_idx_mature = file(params.bwa_idx_mature)
    BWA_MAP_SHORT( fastq_data_s, bwa_idx_mature )
    ch_bam_s = BWA_MAP_SHORT.out.bam_un_short
    SAMTOOLS_FQ_s(ch_bam_s)
    fastq_data_l = channel.fromPath(params.fastq_long).map { file -> tuple(file.baseName, file)}
    bwa_idx_hp = file(params.bwa_idx_hairpin)
    BWA_MAP_LONG( fastq_data_l, bwa_idx_hp )
    ch_bam_l = BWA_MAP_LONG.out.bam_un_long
    SAMTOOLS_FQ_l(ch_bam_l)
}
