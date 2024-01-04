




process create_chunks {
    tag "${fastaFai}"
    publishDir params.outdir+"/parts", mode: 'copy'

    input:
        path fastaFai

    output:
        path '*.bed' , emit: bed

    script:
    if (params.debug){
    """
    echo perl ${baseDir}/scripts/chunk_fai.pl ${fastaFai} ${params.wsize} wpm
    touch wpm1.bed wpm2.bed wpm3.bed
    """
    }else{
    """
    perl ${baseDir}/scripts/chunk_fai.pl ${fastaFai} ${params.wsize} wpm
    """
    }
}


//Run strelka for a multiples samples using a pool!!!!
/*
process STRELKA_POOL{
    tag "$sampleId"
    publishDir "$params.outdir/strelka_pool", mode : "copy"


    conda "bioconda::strelka=2.9.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/strelka:2.9.10--h9ee0642_1' :
        'biocontainers/strelka:2.9.10--h9ee0642_1' }"


    input:
     val(sampleId)
     file(bams)
     file(bais)
     file(part)
    output:
    tuple val("${sampleId}"), path("*variants.vcf.gz")    , emit: vcf
    tuple val("${sampleId}"), path("*variants.vcf.gz.tbi"), emit: vcf_tbi
    tuple val("${sampleId}"), path("genome.*.vcf.gz")      , emit: genome_vcf
    tuple val("${sampleId}"), path("genome.*.vcf.gz.tbi")  , emit: genome_vcf_tbi
    // path ("variants.pass.vcf.gz"), emit: variants


    script:
    def prefix = "${sampleId}"
    def samplesgroup=""
  def bamg=""
    for (b in bams){
        bamg=bamg+" --bam $b"
    }

    if (params.debug == true){
        """
        echo configureStrelkaGermlineWorkflow.py \\
                ${bamg} \\
                --referenceFasta ${params.ref} \\
                --callRegions ${part} \\
                --genome  \\
                --runDir ./strelka_germline

        echo python strelka_germline/runWorkflow.py -m local -j $task.cpus
        touch  genome.${prefix}.vcf.gz genome.${prefix}.vcf.gz.tbi ${prefix}.variants.vcf.gz ${prefix}.variants.vcf.gz.tbi
        """
    }
    else{
        """
     configureStrelkaGermlineWorkflow.py \\
            ${bamg} \\
             --referenceFasta ${params.ref} \\
             --callRegions ${params.brca_reg} \\
             --exome  \\
             --runDir ./strelka_germline
     python strelka_germline/runWorkflow.py -m local -j $task.cpus
     mv strelka_germline/results/variants/genome.*.vcf.gz     .
     mv strelka_germline/results/variants/genome.*.vcf.gz.tbi .
     mv strelka_germline/results/variants/variants.vcf.gz     ${prefix}.variants.vcf.gz
     mv strelka_germline/results/variants/variants.vcf.gz.tbi ${prefix}.variants.vcf.gz.tbi
   """
    }
}

*/

process strelka{
     tag "$bed"
     publishDir "$params.outdir/strelka", mode : "copy"

     conda "bioconda::strelka=2.9.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/strelka:2.9.10--h9ee0642_1' :
        'biocontainers/strelka:2.9.10--h9ee0642_1' }"
    
    input:
      tuple val(sampleId), val(group),file(bams),file(bais),file(bed)
    output:
     tuple val("${part}"), path("*variants.vcf.gz")    , emit: vcf
     tuple val("${part}"), path("*variants.vcf.gz.tbi"), emit: vcf_tbi
     tuple val("${part}"), path("genome.*.vcf.gz")      , emit: genome_vcf
     tuple val("${part}"), path("genome.*.vcf.gz.tbi")  , emit: genome_vcf_tbi
    
    script:
    def part = "all-" + bed.baseName
    def prefix=part
    def samplesgroup=""
    def bamg=""
    for (b in bams){
        bamg=bamg+" --bam $b"
    }
      if (params.debug == true){
        """
        echo configureStrelkaGermlineWorkflow.py \\
                ${bamg} \\
                --referenceFasta ${params.ref} \\
                --callRegions ${bed} \\
                --runDir ./strelka_germline

        echo python strelka_germline/runWorkflow.py -m local -j $task.cpus
        touch  genome.${prefix}.vcf.gz genome.${prefix}.vcf.gz.tbi ${prefix}.variants.vcf.gz ${prefix}.variants.vcf.gz.tbi
        """
    }
    else{
        """
     configureStrelkaGermlineWorkflow.py \\
            ${bamg} \\
             --referenceFasta ${params.ref} \\
             --callRegions ${bed} \\
             --runDir ./strelka_germline
     python strelka_germline/runWorkflow.py -m local -j $task.cpus
     mv strelka_germline/results/variants/genome.*.vcf.gz     .
     mv strelka_germline/results/variants/genome.*.vcf.gz.tbi .
     mv strelka_germline/results/variants/variants.vcf.gz     ${prefix}.variants.vcf.gz
     mv strelka_germline/results/variants/variants.vcf.gz.tbi ${prefix}.variants.vcf.gz.tbi
   """  
    }
}


workflow {
    // TODO do a parameter check
   // PRINT_VERSIONS()
    //we read pairs from regex
 if(params.csv != null){
    //we reads pairs from csv
    read_bams_ch=Channel.fromPath(params.csv) \
        | splitCsv(header:true) \
        | map { row-> tuple(row.sampleId, row.part,  file(row.bam), file(row.bai))}\
    }else{
        println "Error: reads regex or path"
    }

    fai=Channel.fromPath(params.fai)
//    fasta=Channel.fromPath(params.ref)

    all_bamsbais=read_bams_ch.groupTuple(by: 1)

    create_chunks(fai)
    all_parts=create_chunks.out.bed.collect().flatten()
    strparts=all_bamsbais.combine(all_parts)
    strelka(strparts,fasta,fai)

 }   