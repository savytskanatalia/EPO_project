// samtools process > get to .bam; sort; extract unmapped reads

params.sam = "/some/path/*.sam"
params.outdir = "results"


process SAMT_BAM {

    publishDir "${params.outdir}/bam_sorted", mode: 'move', pattern: "*_sorted.bam"

    input:
    tuple val(sample_id), path(sam)


    output:
    path("*_sorted.bam")
    tuple val(sample_id), path("${sample_id}_unmapped.bam"), emit: bam_un


    script:
    """
    samtools view -b $sam | samtools sort  - > "${sample_id}_sorted.bam"
    samtools view  -b -f 4  "${sample_id}_sorted.bam" >  "${sample_id}_unmapped.bam"

    """
}


process SAMTOOLS_FQ {

    publishDir "${params.outdir}/unmapped_fastq", mode: 'move', pattern: "*.fastq.gz"

    input:
    tuple val(sample_id), path(bam)

    output:
    path("*.fastq.gz")

    script:

    """
    samtools fastq "${bam}" > "${sample_id}.fastq.gz"


    """
}


workflow {

    sam_data = channel.fromPath(params.sam).map { file -> tuple(file.baseName, file)}
    SAMT_BAM( sam_data )
    ch_bam = SAMT_BAM.out.bam_un
    SAMTOOLS_FQ(ch_bam)
}
