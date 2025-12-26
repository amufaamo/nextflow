process STAR_INDEX {
    // STAR container definition
    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_8'

    publishDir "${params.star_index_dir}", mode: 'copy', overwrite: false
    cpus params.cpus

    input:
    path fasta
    path gtf

    output:
    path "index", emit: index

    script:
    """
    mkdir -p index
    STAR --runMode genomeGenerate \
        --genomeDir index \
        --genomeFastaFiles ${fasta} \
        --sjdbGTFfile ${gtf} \
        --runThreadN ${task.cpus} \
        --genomeSAindexNbases ${params.star_index_nbases}
    """
}

process STAR_ALIGN {
    tag "$sample_id"
    // ★修正1: バージョンを STAR_INDEX と同じ 2.7.11b に統一しました
    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_8'

    publishDir "${params.outdir}/star/${sample_id}", mode: 'copy'
    cpus params.cpus

    input:
    tuple val(sample_id), path(reads)
    path index
    path gtf

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam
    path "${sample_id}_Log.final.out", emit: log

    script:
    def input_reads = reads.join(' ')
    // ★修正2: パラメータ(例: "30.GB")をバイト数に変換する魔法の呪文です
    def sort_ram = nextflow.util.MemoryUnit.of("${params.star_bam_sort_ram}").toBytes()

    """
    STAR --genomeDir ${index} \
        --readFilesIn ${input_reads} \
        --runThreadN ${task.cpus} \
        --readFilesCommand zcat \
        --outFileNamePrefix ${sample_id}_ \
        --sjdbGTFfile ${gtf} \
        --outSAMtype BAM SortedByCoordinate \
        --limitBAMsortRAM ${sort_ram}
    
    mv ${sample_id}_Aligned.sortedByCoord.out.bam ${sample_id}.bam
    """
}