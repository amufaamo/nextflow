process PIPSEEKER_INDEX {
    executor 'local' 
    cpus 16
    publishDir "${params.outdir}/pipseeker_reference", mode: 'copy'

    input:
    path fasta
    path gtf

    output:
    path "pipseeker_ref_dir", emit: index

    script:
    """
    ${params.pipseeker_bin} buildmapref \
        --fasta ${fasta} \
        --gtf ${gtf} \
        --output-path pipseeker_ref_dir \
        --threads ${task.cpus}
    """
}

process PIPSEEKER_FULL {
    tag "$sample_id"
    cpus 16
    memory '64 GB'
    publishDir "${params.outdir}/pipseeker_output/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path index

    output:
    path "${sample_id}_out", emit: out_dir

    script:
    def r1_name = reads[0].name
    def prefix = r1_name.toString().replaceAll(/_R1.*/, "_")

    """
    mkdir -p ${sample_id}_out
    ${params.pipseeker_bin} full \
        --output-path ${sample_id}_out \
        --fastq ./${prefix} \
        --chemistry ${params.sc_chemistry} \
        --star-index-path ${index} \
        --threads ${task.cpus}
    """
}
