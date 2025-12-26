process BUSCO {
    container 'quay.io/biocontainers/busco:5.5.0--pyhdfd78af_0'
    
    publishDir "${params.outdir}/busco", mode: 'copy'
    cpus params.cpus

    input:
    path fasta
    val lineage // "mammalia_odb10", "vertebrata_odb10" など

    output:
    path "busco_output", emit: dir
    path "busco_output/short_summary.*.txt", emit: summary

    script:
    // ★修正: --offline を削除しました
    """
    busco -i ${fasta} \
          -l ${lineage} \
          -o busco_output \
          -m transcriptome \
          -c ${task.cpus} \
          --download_path ${params.busco_download_path}
    """
}