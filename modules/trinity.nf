process TRINITY {
    // Trinityはメモリを大量に消費するため、必要に応じて調整してください
    container 'quay.io/biocontainers/trinity:2.15.1--h9f5acd7_0'
    
    publishDir "${params.outdir}/trinity", mode: 'copy'
    cpus params.cpus 
    // Trinityには多めのメモリ推奨 (例: 32GB以上)
    memory '32 GB' 

    input:
    path reads_r1 // 全サンプルのRead1のリスト
    path reads_r2 // 全サンプルのRead2のリスト (Singleの場合は空)

    output:
    path "trinity_out_dir/Trinity.fasta", emit: fasta
    path "trinity_out_dir/Trinity.fasta.gene_trans_map", emit: map

    script:
    // ファイルリストをカンマ区切り文字列に変換
    def left_reads  = reads_r1.join(',')
    def right_reads = reads_r2.join(',')
    
    def seqType = "fq"
    def max_memory = "30G" // コンテナ内の制限用

    if (params.single_end) {
        """
        Trinity --seqType ${seqType} \
            --max_memory ${max_memory} \
            --single ${left_reads} \
            --CPU ${task.cpus} \
            --output trinity_out_dir \
            --full_cleanup
        """
    } else {
        """
        Trinity --seqType ${seqType} \
            --max_memory ${max_memory} \
            --left ${left_reads} \
            --right ${right_reads} \
            --CPU ${task.cpus} \
            --output trinity_out_dir \
            --full_cleanup
        """
    }
}
