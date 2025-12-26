process TRINITY {
    container 'trinityrnaseq/trinityrnaseq:2.14.0'
    
    publishDir "${params.outdir}/trinity", mode: 'copy', overwrite: true
    
    cpus params.cpus 
    memory '200 GB' 

    input:
    path reads_r1 
    path reads_r2 

    output:
    // リネーム後の標準的なファイル名で出力します
    path "Trinity.fasta", emit: fasta
    path "Trinity.fasta.gene_trans_map", emit: map

    script:
    def left_reads  = reads_r1.join(',')
    def right_reads = reads_r2.join(',')
    
    def seqType = "fq"
    def max_memory = "180G" 

    // コマンド実行後に mv コマンドでファイル名を整形します
    if (params.single_end) {
        """
        Trinity --seqType ${seqType} \
            --max_memory ${max_memory} \
            --single ${left_reads} \
            --CPU ${task.cpus} \
            --output trinity_out_dir \
            --full_cleanup
        
        mv trinity_out_dir.Trinity.fasta Trinity.fasta
        mv trinity_out_dir.Trinity.fasta.gene_trans_map Trinity.fasta.gene_trans_map
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

        mv trinity_out_dir.Trinity.fasta Trinity.fasta
        mv trinity_out_dir.Trinity.fasta.gene_trans_map Trinity.fasta.gene_trans_map
        """
    }
}