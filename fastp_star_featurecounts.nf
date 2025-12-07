#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * ================================================================================================
 * パイプラインのパラメータ定義
 * ================================================================================================
 */
// --- 入力/出力ファイルのパス ---
params.samplesheet = './samples.csv'                  // サンプルシートのパス
params.fastq_dir = "${System.getProperty("user.home")}/fastq" // FASTQファイルが保存されているディレクトリ
params.ref_fasta = "path/to/your/genome.fna"           // リファレンスゲノムFASTA
params.ref_gtf = "path/to/your/genome.gtf"             // リファレンスアノテーションGTF
params.adapter_fasta = "${System.getProperty("user.home")}/fasta/adapter.fasta" // fastpで使うアダプターFASTA
params.star_index = "./star_index"                     // STARのインデックスを保存するディレクトリ
params.results_folder_name = "results"                 // 結果を保存するフォルダの名前
params.outdir = "txt/${params.results_folder_name}"    // 結果を出力する親ディレクトリ
params.cpus = 8                                        // 使用するCPUスレッド数

// --- 各ツールの実行パス ---
params.fastp_path = "/home/masakazu/miniconda3/envs/fastp/bin/fastp"
params.star_path = "/home/masakazu/miniconda3/envs/star/bin/STAR"
params.featurecounts_path = "/home/masakazu/miniconda3/envs/featurecounts/bin/featureCounts"


// パイプライン開始時にパラメータをログ表示
log.info """
         R N A - S E Q   P I P E L I N E (Simplified ver.)
         ================================================
         SampleSheet   : ${params.samplesheet}
         FASTQ Dir     : ${params.fastq_dir}
         Reference FASTA: ${params.ref_fasta}
         Reference GTF : ${params.ref_gtf}
         Adapter FASTA : ${params.adapter_fasta}
         Output dir    : ${params.outdir}
         ================================================
         """
         .stripIndent()

/*
 * ================================================================================================
 * ワークフローの本体
 * ================================================================================================
 */
workflow {
    Channel
        .fromPath( params.samplesheet )
        .splitCsv( header:true )
        .map { row ->
            def sample_id = row.sample
            def reads = [
                file("${params.fastq_dir}/${row.fastq_1}"),
                file("${params.fastq_dir}/${row.fastq_2}")
            ]
            tuple( sample_id, reads )
        }
        .ifEmpty { error "Sample sheet '${params.samplesheet}' is empty, not found, or contains file paths that could not be found in '${params.fastq_dir}'" }
        .set { ch_reads }

    // --- シンプルになったワークフロー ---
    FASTP( ch_reads )
    STAR_INDEX( params.ref_fasta, params.ref_gtf )
    STAR_ALIGN( FASTP.out.reads, STAR_INDEX.out.index )
    FEATURECOUNTS( STAR_ALIGN.out.bam, params.ref_gtf )
}

/*
 * ================================================================================================
 * プロセスの定義
 * ================================================================================================
 */

process FASTP {
    publishDir "${params.outdir}/fastp/${sample_id}", mode: 'copy'
    input:
    tuple val(sample_id), path(reads)
    output:
    tuple val(sample_id), path("${sample_id}_Read{1,2}_trimmed.fastq.gz"), emit: reads
    path "${sample_id}.fastp.html", emit: html
    path "${sample_id}.fastp.json", emit: json
    script:
    def read1_out = "${sample_id}_Read1_trimmed.fastq.gz"
    def read2_out = "${sample_id}_Read2_trimmed.fastq.gz"
    """
    ${params.fastp_path} \\
        -i ${reads[0]} \\
        -I ${reads[1]} \\
        -o ${read1_out} \\
        -O ${read2_out} \\
        -h ${sample_id}.fastp.html \\
        -j ${sample_id}.fastp.json \\
        --thread ${task.cpus} \\
        --trim_poly_g \\
        --trim_poly_x \\
        --adapter_fasta ${params.adapter_fasta}
    """
}

process STAR_INDEX {
    publishDir "${params.star_index}", mode: 'copy', overwrite: false
    input:
    path fasta
    path gtf
    output:
    path "index", emit: index
    script:
    """
    ${params.star_path} \\
        --runMode genomeGenerate \\
        --genomeDir index \\
        --genomeFastaFiles ${fasta} \\
        --sjdbGTFfile ${gtf} \\
        --runThreadN ${task.cpus}
    """
}

process STAR_ALIGN {
    publishDir "${params.outdir}/star/${sample_id}", mode: 'copy'
    input:
    tuple val(sample_id), path(reads)
    path index
    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam
    path "${sample_id}_Log.final.out", emit: log
    script:
    """
    ${params.star_path} \\
        --genomeDir ${index} \\
        --readFilesIn ${reads[0]} ${reads[1]} \\
        --runThreadN ${task.cpus} \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${sample_id}_ \\
        --outSAMtype BAM SortedByCoordinate
    
    mv ${sample_id}_Aligned.sortedByCoord.out.bam ${sample_id}.bam
    """
}

process FEATURECOUNTS {
    publishDir "${params.outdir}/featurecounts", mode: 'copy'
    input:
    tuple val(sample_id), path(bam)
    path gtf
    output:
    path "${sample_id}_featurecounts.txt"
    path "${sample_id}_featurecounts.txt.summary"
    script:
    """
    ${params.featurecounts_path} \\
        -p \\
        --countReadPairs \\
        -a ${gtf} \\
        -o ${sample_id}_featurecounts.txt \\
        -T ${task.cpus} \\
        ${bam}
    """
}