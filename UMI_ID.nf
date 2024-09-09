//  reads three files with similar IDs - R1,R2 and index >>> to extract UMI from index and put it in read IDs
// "/path/*{index,R1,R2}_001.fastq.gz" 
// from path or from pairs? channel.fromFilePairs
// as fastp tries to work in pairs to utilize it I process R1+index and R2+index as separate instances
// need to use a single core only
// possible improvement: implement pecheck?! just in case to make sure the pairs are in synch
// --disable_adapter_trimming --disable_quality_filtering 
// --trim_poly_g because we had nova_seq 


// Example Usage: 
// nextflow run UMI_ID.nf -c generic.config --my_files "/mydir/*{index,R1,R2}_001.fastq.gz" -with-docker savytskanatalia/quantefication2 --outdir "results"



params.my_files = "/path/*{index,R1,R2}_001.fastq.gz" 
params.outdir = "results"


process FASTP_UMI {
    publishDir 'results'
    publishDir "${params.outdir}/trimmed", mode: 'move', pattern: "*_PROC.fastq.gz"
    publishDir "${params.outdir}/log", mode: 'move', pattern: "*_fastp.json"
    publishDir "${params.outdir}/log", mode: 'move', pattern: "*_fastp.html"
    publishDir "${params.outdir}/failed", mode: 'move', pattern: "*_failed.fastq.gz"
    input:
    tuple val(sampleId), file(reads)
    output:
    path("*_R{1,2}_PROC.fastq.gz")
    path("*_fastp.json")
    path("*_fastp.html")
    path("*_failed.fastq.gz")
    script:
    """
    echo HI THERE ${reads[0]} ${reads[1]} ${reads[2]} $sampleId
    awk -v FS="\t" -v OFS="\t" 'NR==FNR {split(\$1, id, " "); umi[id[1]]=\$2;  next;} {split(\$1, id, " "); \$1=id[1]":"umi[id[1]]" "id[2]; print \$0}'  <(zcat ${reads[2]}|paste - - - -) <(zcat ${
reads[0]} |paste - - - -)|tr "\t" "\n"|bgzip -c > ${sampleId}_R1_UMI.fastq.gz
    awk -v FS="\t" -v OFS="\t" 'NR==FNR {split(\$1, id, " "); umi[id[1]]=\$2;  next;} {split(\$1, id, " "); \$1=id[1]":"umi[id[1]]" "id[2]; print \$0}'  <(zcat ${reads[2]}|paste - - - -) <(zcat ${
reads[1]} |paste - - - -)|tr "\t" "\n"|bgzip -c > ${sampleId}_R2_UMI.fastq.gz
    fastp -i ${sampleId}_R1_UMI.fastq.gz -I ${sampleId}_R2_UMI.fastq.gz -o ${sampleId}_R1_PROC.fastq.gz -O ${sampleId}_R2_PROC.fastq.gz --detect_adapter_for_pe -w 1 -p -j ${sampleId}_fastp.json -h
  ${sampleId}_fastp.html  --failed_out ${sampleId}_failed.fastq.gz --trim_poly_g --disable_quality_filtering 

    """
}





workflow {
  fastq_files = Channel.fromFilePairs(params.my_files, size: -1)
  FASTP_UMI(fastq_files)
}

