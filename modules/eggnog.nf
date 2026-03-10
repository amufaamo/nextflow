process EGGNOG_MAPPER {
    container 'quay.io/biocontainers/eggnog-mapper:2.1.9--pyhdfd78af_0'
    publishDir "${params.outdir}/annotation/eggnog", mode: 'copy'
    cpus params.cpus

    input:
    path pep_fasta
    path eggnog_db_dir

    output:
    path "*.emapper.annotations", emit: annotations
    path "*.emapper.hits", emit: hits

    script:
    """
    emapper.py -i ${pep_fasta} --itype proteins -o eggnog_out -m diamond --data_dir ${eggnog_db_dir} --cpu ${task.cpus}
    """
}
