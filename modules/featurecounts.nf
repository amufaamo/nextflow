process FEATURECOUNTS {
    tag "$sample_id"
    // Subread container definition
    container 'quay.io/biocontainers/subread:2.0.1--h5bf99c6_1'

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
    
    // もしパラメータが空なら 'gene_id' を代入する変数を作る
    def group_feature = params.fc_group_features ?: 'gene_id'

    """
    featureCounts \
        ${paired_opts} \
        -t exon \
        -g ${group_feature} \
        -a ${gtf} \
        -o ${sample_id}_featurecounts.txt \
        -T ${task.cpus} \
        ${bam}
    """
}