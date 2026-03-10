process DESEQ2 {
    container 'quay.io/biocontainers/mulled-v2-886915b4ba251912df665aef34685ff87e743a18:595015b36fa2c129e7019f864ea6526ea02c525f-0'
    publishDir "${params.outdir}/deseq2", mode: 'copy'

    input:
    path counts_txt

    output:
    path "deseq2_results.csv", emit: results

    script:
    """
    #!/usr/bin/env Rscript
    library(DESeq2)
    
    # Corsetカウントファイルの読み込み
    counts <- read.table("${counts_txt}", header=TRUE, row.names=1)
    
    # 簡易的なメタデータの自動付与 (前半グループ vs 後半グループ)
    # 実運用では外部の metadata.csv を読み込む形に拡張できます
    samples <- colnames(counts)
    group <- factor(rep(c("Condition_A", "Condition_B"), length.out=length(samples)))
    
    colData <- data.frame(row.names=samples, condition=group)
    dds <- DESeqDataSetFromMatrix(countData=counts, colData=colData, design=~condition)
    dds <- DESeq(dds)
    res <- results(dds)
    
    write.csv(as.data.frame(res), file="deseq2_results.csv")
    """
}
