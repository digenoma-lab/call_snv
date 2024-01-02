



process create_chunks {
    tag "${fastaFai}"
    publishDir params.outdir+"/parts", mode: 'copy'

    input:
        file(fastaFai)

    output:
        path '*.bed' , emit: bed

    script:
    if (params.debug){
    """
    echo perl ${baseDir}/scripts/split-faix.pl ${fastaFai} ${params.nsplit}
    touch pp1.bed pp2.bed pp3.bed
    """
    }else{
    """
    perl ${baseDir}/scripts/split-faix.pl ${fastaFai} ${params.nsplit}
    """
    }
}


workflow {




 }   