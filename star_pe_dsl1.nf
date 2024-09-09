#!/usr/bin/env nextflow
nextflow.enable.dsl=1
/*
================================================================================
                           STAR-mapping-PE
================================================================================
  Nextflow pipeline for RNA-seq mapping with STAR v2.7.5
--------------------------------------------------------------------------------
 @Repository
 https://github.com/savytskanatalia/
 @Author
 Savytska Natalia, 2024. Genome Biology of Neurodegenerative Diseases, DZNE TU
--------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run star_pe_dsl1.nf -c star_pe_dsl1.config --reads "/absolute/path/to/files/*_R{1,2}.fastq.gz" ... -with-docker savytskanatalia/quantefication2 
    Mandatory arguments:
      --reads [str]                  Path to input data (must be surrounded with quotes). Default: "input/*_R{1,2}*.fastq.gz" 
      --outdir [str]                 Output directory path. Default: output/      
    References                        If not specified in the configuration file or you wish to overwrite any of the references.
      --star_index [str]             Path to STAR-Index reference
      --staropt [str]                Additional STAR flags for processing. If any additional STAR parameters need to be included, they can be with this flag (e.g. '--alignIntronMin 20 ...').
    Options:
      --threads [int]                Number of threads to use. Default: 4
      --read_length [int]            Length of the reads. Default: 100
      --mism [int]                   outFilterMismatchNmax in STAR  alignment. Default: 999
      --misml [int]                  outFilterMismatchNoverLmax in STAR  alignment. Default: 0.1
      --sorted [str]                 Option for STAR and TElocal. If .bam files are sorted. Default: 1 (=unsorted, no additional parameter for TElocal run). Other valid option: 2 (sorted)
    """.stripIndent()
}


sortt=params.sorted
fastq_files=Channel.fromFilePairs(params.reads)
out=params.outdir
read_length=Channel.value(params.read_length)
threads=Channel.value(params.threads)
mismatch_n=Channel.value(params.mism)
mismatch_nl=Channel.value(params.misml)
outdir=params.outdir
staropt=params.staropt



process star_mapping {
        container='savytskanatalia/quantefication2'
        echo true
        publishDir "${params.outdir}/mapped", mode: 'copy'
        input:
        set sampleId, file(reads) from fastq_files
        path x from Channel.value(params.star_index)
        val y from threads
        val a from mismatch_n
        val b from mismatch_nl
        val s from staropt
        output:
        file "*.bam"
        file "*"
        script:
        if( sortt == 1)
            """
            source activate telocal
            pwd
            echo HI THERE ${reads[0]} $sampleId
            STAR \\
                --outSAMstrandField intronMotif \\
                --outSAMunmapped Within \\
                --runThreadN $y \\
                --genomeDir $x \\
                --readFilesIn <(zcat ${reads[0]}) <(zcat ${reads[1]})  \\
                --outFilterMultimapNmax 100  \\
                --winAnchorMultimapNmax 100 \\
                --outSAMtype BAM Unsorted \\
                --outFilterMismatchNmax $a \\
                --outFilterMismatchNoverLmax $b \\
                --outFileNamePrefix  ${sampleId}_MM_100 $s


            """
        else if(sortt == 2)
            """
            source activate telocal
            pwd
            echo HI THERE ${reads[0]} $sampleId
            STAR \\
                --outSAMstrandField intronMotif \\
                --outSAMunmapped Within \\
                --runThreadN $y \\
                --genomeDir $x \\
                --readFilesIn <(zcat ${reads[0]}) <(zcat ${reads[1]}) \\
                --outFilterMultimapNmax 100  \\
                --winAnchorMultimapNmax 100 \\
                --outSAMtype BAM SortedByCoordinate \\
                --outFilterMismatchNmax $a \\
                --outFilterMismatchNoverLmax $b \\
                --outFileNamePrefix  ${sampleId}_MM_100 $s


            """
        else
            error "Invalid alignment mode: ${sortt}"
}
