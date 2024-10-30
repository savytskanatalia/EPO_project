                                                                                                             
// using bbduk to split single-end reads by length (31) into separate files
// part of the miRNA/small RNA analysis



params.fastq = "/some/path/*.fastq.gz"
params.outdir = "results"


process BBDUK_SPLIT {
    publishDir "${params.outdir}/longerReads", mode: 'move', pattern: "*_longReads.fastq.gz"
    publishDir "${params.outdir}/below31bpReads", mode: 'move', pattern: "*_shortReads.fastq.gz"

    input:
    path fq

    output:
    path("*_longReads.fastq.gz")
    path("*_shortReads.fastq.gz")

    script:
    """
    bbduk.sh in=$fq out=${fq.baseName}_longReads.fastq minlen=32 outm=${fq.baseName}_shortReads.fastq
    gzip ${fq.baseName}_longReads.fastq
    gzip ${fq.baseName}_shortReads.fastq

    """
}



workflow {
  fq_files = Channel.fromPath(params.fastq)
  BBDUK_SPLIT(fq_files)
}
