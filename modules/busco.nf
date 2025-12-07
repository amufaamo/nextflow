process BUSCO {
    container 'quay.io/biocontainers/busco:5.5.0--pyhdfd78af_0'
    
    publishDir "${params.outdir}/busco", mode: 'copy'
    cpus params.cpus

    input:
    path fasta
    val lineage // "mammalia_odb10", "eukaryota_odb10" など

    output:
    path "busco_output", emit: dir
    path "busco_output/short_summary.*.txt", emit: summary

    script:
    """
    busco -i ${fasta} \
          -l ${lineage} \
          -o busco_output \
          -m transcriptome \
          -c ${task.cpus} \
          --offline \
          --download_path ${params.busco_download_path} 
    
    # 注: 初回実行時にデータセットがない場合、--offlineを外すか
    # 事前にダウンロードしておく必要があります。
    # 今回はシンプルにするため自動ダウンロードを許可する設定にします:
    # busco -i ${fasta} -l ${lineage} -o busco_output -m transcriptome -c ${task.cpus}
    """
}
