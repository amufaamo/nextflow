process CORSET {
    container 'quay.io/biocontainers/corset:1.09--h9f5acd7_3'
    publishDir "${params.outdir}/corset", mode: 'copy'
    cpus params.cpus

    input:
    path eq_classes // salmonのeq_classes.txt のリスト

    output:
    path "counts.txt", emit: counts
    path "clusters.txt", emit: clusters

    script:
    """
    # 複数渡された eq_classes.txt を解釈して実行
    corset -p ${task.cpus} -i salmon ${eq_classes.join(' ')}
    """
}
