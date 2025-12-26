process FEATURECOUNTS {
    tag "$sample_id"
    // Subread container definition (v2.1.1 なので --countReadPairs が使えます)
    container 'quay.io/biocontainers/subread:2.1.1--h577a1d6_0'

    publishDir "${params.outdir}/featurecounts", mode: 'copy'
    cpus params.cpus

    input:
    tuple val(sample_id), path(bam)
    path gtf

    output:
    path "${sample_id}_featurecounts.txt"
    path "${sample_id}_featurecounts.txt.summary"

    script:
    def paired_opts = params.single_end ? "" : "-p --countReadPairs"
    
    // パラメータのデフォルト値処理
    def group_feature = params.fc_group_features ?: 'gene_id'
    
    // ★ここが修正ポイント！ configで定義したパラメータを使います
    def feature_type = params.gtf_feature_type ?: 'exon'

    """
    featureCounts \
        ${paired_opts} \
        -t ${feature_type} \
        -g ${group_feature} \
        -a ${gtf} \
        -o ${sample_id}_featurecounts.txt \
        -T ${task.cpus} \
        ${bam}
    """
}