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
    path "${sample_id}_quant/aux_info/eq_classes.txt", emit: eq_classes

    script:
    def is_single = reads instanceof Path || reads.size() == 1
    if (is_single) {
        def in_read1 = reads instanceof Path ? reads : reads[0]
        """
        salmon quant -i ${index} -l A \\
            -r ${in_read1} \\
            -p ${task.cpus} \\
            --validateMappings \\
            --dumpEq \\
            -o ${sample_id}_quant
        """
    } else {
        def in_read1 = reads[0]
        def in_read2 = reads[1]
        """
        salmon quant -i ${index} -l A \\
            -1 ${in_read1} -2 ${in_read2} \\
            -p ${task.cpus} \\
            --validateMappings \\
            --dumpEq \\
            -o ${sample_id}_quant
        """
    }
}
