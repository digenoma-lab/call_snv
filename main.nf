




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

//Multisample strelka SNP calling
process strelka {
     tag "all"
     publishDir "$params.outdir/strelka", mode : "copy"

    conda "bioconda::strelka=2.9.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/strelka:2.9.10--h9ee0642_1' :
        'biocontainers/strelka:2.9.10--h9ee0642_1' }" 
    
    input:
      tuple val(sampleId), val(group),file(bams),file(bais)
    output:
     tuple val("${part}"), path("*variants.vcf.gz")    , emit: vcf
     tuple val("${part}"), path("*variants.vcf.gz.tbi"), emit: vcf_tbi
     tuple val("${part}"), path("genome.*.vcf.gz")      , emit: genome_vcf
     tuple val("${part}"), path("genome.*.vcf.gz.tbi")  , emit: genome_vcf_tbi
    
    script:
    def part = "all-"
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
                --callRegions ${params.bed} \\
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
             --callRegions ${params.bed} \\
             --runDir ./strelka_germline
     python strelka_germline/runWorkflow.py -m local -j $task.cpus
     mv strelka_germline/results/variants/genome.*.vcf.gz     .
     mv strelka_germline/results/variants/genome.*.vcf.gz.tbi .
     mv strelka_germline/results/variants/variants.vcf.gz     ${prefix}.variants.vcf.gz
     mv strelka_germline/results/variants/variants.vcf.gz.tbi ${prefix}.variants.vcf.gz.tbi
   """  
    }
}



//Multisample manta SV calling
process manta {
     tag "all"
     publishDir "$params.outdir/manta", mode : "copy" 
    
    input:
      tuple val(sampleId), val(group),file(bams),file(bais)
    output:    
  tuple val(part), path("manta/Manta_${part}.candidateSmallIndels.vcf.gz")    , emit: manta_small_indels
  tuple val(part), path("manta/Manta_${part}.candidateSmallIndels.vcf.gz.tbi")  , emit: manta_small_indels_tbi
  tuple val(part), path("manta/Manta_${part}.candidateSV.vcf.gz")         , emit: manta_candidate
  tuple val(part), path("manta/Manta_${part}.candidateSV.vcf.gz.tbi")       , emit: manta_candidate_tbi
  tuple val(part), path("manta/Manta_${part}.diploidSV.vcf.gz")         , emit: manta_diploid
  tuple val(part), path("manta/Manta_${part}.diploidSV.vcf.gz.tbi")       , emit: manta_diploid_tbi
  tuple val(part), path("manta/Manta_${part}.diploidSV_converted.vcf.gz")     , emit: manta_diploid_convert
  tuple val(part), path("manta/Manta_${part}.diploidSV_converted.vcf.gz.tbi")   , emit: manta_diploid_convert_tbi


    script:
    def part = "all"
    def prefix=part
    def samplesgroup=""
    def bamg=""
    for (b in bams){
        bamg=bamg+" --bam $b"
    }
      if (params.debug == true){
        """
        echo configManta.py  \\
                ${bamg} \\
                --referenceFasta ${params.ref} \\
                --callRegions ${params.bed} \\
                --runDir ./manta

        echo python strelka_germline/runWorkflow.py -m local -j $task.cpus
        mkdir -p manta
        touch manta/Manta_${prefix}.candidateSmallIndels.vcf.gz  manta/Manta_${prefix}.candidateSmallIndels.vcf.gz.tbi
        touch manta/Manta_${prefix}.candidateSV.vcf.gz manta/Manta_${prefix}.candidateSV.vcf.gz.tbi
        touch manta/Manta_${prefix}.diploidSV.vcf.gz manta/Manta_${prefix}.diploidSV.vcf.gz.tbi 
        touch manta/Manta_${sampleID}.diploidSV_converted.vcf manta/Manta_${part}.diploidSV_converted.vcf.gz.tbi
        """
    }
    else{
        """
     configManta.py  \\
            ${bamg} \\
             --referenceFasta ${params.ref} \\
             --callRegions ${params.bed} \\
             --runDir ./manta
    manta/runWorkflow.py -m local -j ${task.cpus}     
    # clean up outputs
  mv manta/results/variants/candidateSmallIndels.vcf.gz \
    manta/Manta_${prefix}.candidateSmallIndels.vcf.gz
  mv manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
    manta/Manta_${prefix}.candidateSmallIndels.vcf.gz.tbi
  mv manta/results/variants/candidateSV.vcf.gz \
    manta/Manta_${prefix}.candidateSV.vcf.gz
  mv manta/results/variants/candidateSV.vcf.gz.tbi \
    manta/Manta_${prefix}.candidateSV.vcf.gz.tbi
  mv manta/results/variants/diploidSV.vcf.gz \
    manta/Manta_${prefix}.diploidSV.vcf.gz
  mv manta/results/variants/diploidSV.vcf.gz.tbi \
    manta/Manta_${prefix}.diploidSV.vcf.gz.tbi  

# convert multiline inversion BNDs from manta vcf to single line
  convertInversion.py \$(which samtools) ${params.ref} \
    manta/Manta_${sampleID}.diploidSV.vcf.gz \
    > manta/Manta_${sampleID}.diploidSV_converted.vcf

  # zip and index converted vcf
  bgzip manta/Manta_${sampleID}.diploidSV_converted.vcf
  tabix manta/Manta_${sampleID}.diploidSV_converted.vcf.gz

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

    all_bamsbais=read_bams_ch.groupTuple(by: 1)

    /*create_chunks(fai)
    all_beds=create_chunks.out.bed.collect().flatten()
    all_beds=all_beds.map { [it.baseName,
          it] }

    all_tbi=create_chunks.out.tbi.collect().flatten()
    all_tbi=all_tbi.map { [it.baseName.replaceAll(".gz",""), it]}

    all_beds=all_beds.combine(all_tbi,by: 0)
    strparts=all_bamsbais.combine(all_beds)

    strparts.view()
    */

    if (params.runstrelka){
    strelka(all_bamsbais)
    }

    if(params.runmanta){
      manta(all_bamsbais)
    }

 }   