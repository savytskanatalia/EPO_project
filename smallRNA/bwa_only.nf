// only bwa mapping
// -with-docker 3a962903b45b
params.fastq_short = "/some/path/*.fastq.gz"
params.fastq_long = "/some/path/*.fastq.gz"
params.bwa_idx_mature = "data/genome.fa.{,amb,ann,bwt,pac,sa}"
params.bwa_idx_hairpin = "data/genome.fa.{,amb,ann,bwt,pac,sa}"
params.outdir = "results"


process BWA_MAP_SHORT {
    publishDir "${params.outdir}/mapped/below31bpReads", mode: 'move', pattern: "*.sam"

    input:
    tuple val(sample_id), path(fastq)
    path bwa_idx_mature

    output:
    path("*.sam")


    script:
    def idxbase = bwa_idx_mature[0].baseName
    """
    bwa mem -t 3 ${idxbase} ${fastq} > ${sample_id}.sam

    """
}
process BWA_MAP_LONG {
    publishDir "${params.outdir}/mapped/longerReads", mode: 'move', pattern: "*.sam"

    input:
    tuple val(sample_id), path(fastq)
    path bwa_idx_mature

    output:
    path("*.sam")

    script:
    def idxbase = bwa_idx_mature[0].baseName
    """
    bwa mem -t 3 ${idxbase} ${fastq} > ${sample_id}.sam
    """
}


workflow {

    fastq_data_s = channel.fromPath(params.fastq_short).map { file -> tuple(file.baseName, file)}
    bwa_idx_mature = file(params.bwa_idx_mature)
    BWA_MAP_SHORT(fastq_data_s, bwa_idx_mature )

    fastq_data_l = channel.fromPath(params.fastq_long).map { file -> tuple(file.baseName, file)}
    bwa_idx_hp = file(params.bwa_idx_hairpin)
    BWA_MAP_LONG( fastq_data_l, bwa_idx_hp )

}
