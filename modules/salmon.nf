process SALMON_INDEX {
    container 'quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0'
    
    publishDir "${params.outdir}/salmon/index", mode: 'copy'
    cpus params.cpus

    input:
    path transcriptome_fasta

    output:
    path "salmon_index", emit: index

    script:
    """
    salmon index -t ${transcriptome_fasta} -i salmon_index -p ${task.cpus}
    """
}

process SALMON_QUANT {
    tag "$sample_id"
    container 'quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0'
    
    publishDir "${params.outdir}/salmon/${sample_id}", mode: 'copy'
    cpus params.cpus

    input:
    tuple val(sample_id), path(reads)
    path index

    output:
    path "${sample_id}_quant", emit: quant_dir

    script:
    if (params.single_end) {
        """
        salmon quant -i ${index} -l A \
            -r ${reads[0]} \
            -p ${task.cpus} \
            --validateMappings \
            -o ${sample_id}_quant
        """
    } else {
        """
        salmon quant -i ${index} -l A \
            -1 ${reads[0]} -2 ${reads[1]} \
            -p ${task.cpus} \
            --validateMappings \
            -o ${sample_id}_quant
        """
    }
}
