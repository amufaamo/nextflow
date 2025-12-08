process TRINITY {
    // ★修正: 確実に存在する安定版(2.14.0)に変更
    container 'trinityrnaseq/trinityrnaseq:2.14.0'
    
    publishDir "${params.outdir}/trinity", mode: 'copy', overwrite: true
    
    // マシンスペックに合わせてメモリ最大化
    cpus params.cpus 
    memory '200 GB' 

    input:
    path reads_r1 
    path reads_r2 

    output:
    path "trinity_out_dir/Trinity.fasta", emit: fasta
    path "trinity_out_dir/Trinity.fasta.gene_trans_map", emit: map

    script:
    def left_reads  = reads_r1.join(',')
    def right_reads = reads_r2.join(',')
    
    def seqType = "fq"
    def max_memory = "180G" 

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