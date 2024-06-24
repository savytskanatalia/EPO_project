//  reads three files with similar IDs - R1,R2 and index >>> to extract UMI from index and put it in read IDs
// "/path/*{index,R1,R2}_001.fastq.gz" 
// from path or from pairs? channel.fromFilePairs
// as fastp tries to work in pairs to utilize it I process R1+index and R2+index as separate instances
// need to use a single core only
// implement pecheck?! just in case to make sure the pairs are in synch
// make it two runs - first will extract UMIs, second take care of adapters in proper pairs with extracted UMIs; we do NOT apply control filtering
// --disable_adapter_trimming --disable_quality_filtering 
// https://github.com/OpenGene/fastp/issues/23

params.my_files = "/path/*{index,R1,R2}_001.fastq.gz" 



process FASTP_UMI {
    input:
    tuple val(sampleId), file(reads)
    output:
    path("*.bam"), emit: bam


    script:
    """
    echo HI THERE ${reads[0]} $sampleId
    fastp -i ${reads[0]} -I ${reads[1]} -o ${sampleId}_R1_PROC.fastq.gz -O ${sampleId}_R2_PROC.fastq.gz --detect_adapter_for_pe -w 1 -p -j ${sampleId}_fastp.json -h  ${sampleId}_fastp.html  --failed_out ${sampleId}_failed.fastq.gz

    fastp -i ${reads[0]} -I ${reads[2]} -o ${sampleId}_R1_PROC.fastq.gz -O ${sampleId}_trash1.out.fq --umi --umi_loc=read2 --umi_len=12 -Q -A -L -w 1
    fastp -i ${reads[1]} -I ${reads[2]} -o ${sampleId}_R2_PROC.fastq.gz -O ${sampleId}_trash2.out.fq --umi --umi_loc=read2 --umi_len=12 -Q -A -L -w 1

    """
}




workflow {
  sam_files = Channel.fromPath(params.sam)
  SAMTOOLS_VIEW_BAM(sam_files)
  ch_bam = SAMTOOLS_VIEW_BAM.out.bam
  SAMTOOLS_SORT_INDEX(ch_bam)
}
