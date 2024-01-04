




process create_chunks {
    tag "${fastaFai}"
    publishDir params.outdir+"/parts", mode: 'copy'

    input:
        path fastaFai

    output:
       path '*.bed.gz' , emit: bed
       path '*.bed.gz.tbi' , emit: tbi 

    script: 
    if (params.debug){
    """
    echo perl ${baseDir}/scripts/chunk_fai.pl ${fastaFai} ${params.wsize} wpm
    touch wpm1.bed.gz wpm2.bed.gz wpm3.bed.gz
    touch wpm1.bed.gz.tbi wpm2.bed.gz.tbi wpm3.bed.gz.tbi
    """
    }else{
    """
    perl ${baseDir}/scripts/chunk_fai.pl ${fastaFai} ${params.wsize} wpm
    ls *.bed | xargs -n 1 bgzip
    ls *.gz | xargs -n 1 tabix 
    """
    }
}


process strelka {
     tag "$bed"
     publishDir "$params.outdir/strelka", mode : "copy"

    conda "bioconda::strelka=2.9.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/strelka:2.9.10--h9ee0642_1' :
        'biocontainers/strelka:2.9.10--h9ee0642_1' }" 
    
    input:
      tuple val(sampleId), val(group),file(bams),file(bais),val(bed),file(b),file(t)
    output:
     tuple val("${part}"), path("*variants.vcf.gz")    , emit: vcf
     tuple val("${part}"), path("*variants.vcf.gz.tbi"), emit: vcf_tbi
     tuple val("${part}"), path("genome.*.vcf.gz")      , emit: genome_vcf
     tuple val("${part}"), path("genome.*.vcf.gz.tbi")  , emit: genome_vcf_tbi
    
    script:
    def part = "all-" + bed.replaceAll(".bed","")
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
                --callRegions ${b} \\
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
             --callRegions ${b} \\
             --runDir ./strelka_germline
     python strelka_germline/runWorkflow.py -m local -j $task.cpus
     mv strelka_germline/results/variants/genome.*.vcf.gz     .
     mv strelka_germline/results/variants/genome.*.vcf.gz.tbi .
     mv strelka_germline/results/variants/variants.vcf.gz     ${prefix}.variants.vcf.gz
     mv strelka_germline/results/variants/variants.vcf.gz.tbi ${prefix}.variants.vcf.gz.tbi
   """  
    }
}


workflow{
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
    all_beds=create_chunks.out.bed.collect().flatten()
    all_beds=all_beds.map { [it.baseName,
          it] }

    all_tbi=create_chunks.out.tbi.collect().flatten()
    all_tbi=all_tbi.map { [it.baseName.replaceAll(".gz",""), it]}

    all_beds=all_beds.combine(all_tbi,by: 0)
    strparts=all_bamsbais.combine(all_beds)
    strparts.view()
    strelka(strparts)

 }   