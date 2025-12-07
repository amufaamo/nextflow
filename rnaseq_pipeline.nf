#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * ================================================================================================
 * 🌟 Universal RNA-Seq Pipeline 🌟
 * バージョン表示・ncRNA対応・Index自動作成ロジック完全版
 * ================================================================================================
 */

// ★★★ ここでバージョンを管理！ ★★★
def pipeline_version = "1.2.0"

// --- 必須パラメータ ---
// GTFは定量(featureCounts)で必ず必要です。ncRNA検出の要です。
params.ref_gtf     = null 

// --- インデックス関連（どちらか片方が必須）---
// 既存のインデックスを使うなら --star_index
// 新しく作るなら --ref_fasta を指定（--star_index を指定しなければ自動で作ります）
params.ref_fasta   = null
params.star_index  = null

// --- オプションパラメータ ---
params.samplesheet = './samples.csv'
params.fastq_dir   = "${System.getProperty("user.home")}/fastq"
params.outdir      = "results"
params.cpus        = 8

// --- ツールのパス設定 ---
params.fastp_path         = 'fastp'
params.star_path          = 'STAR'
params.featurecounts_path = 'featureCounts'
params.adapter_fasta      = "${System.getProperty("user.home")}/fasta/adapter.fasta" 
params.star_index_dir     = "./star_index_new" // 新規作成時の保存先名

// --- ツールバージョン取得関数 ---
def get_version(cmd_path, version_flag) {
    try {
        def proc = [cmd_path, version_flag].execute()
        proc.waitFor()
        def out = proc.in.text.trim()
        def err = proc.err.text.trim()
        def ver_str = (out + err).readLines().find { it.trim() != "" } ?: "Unknown"
        return ver_str
    } catch (Exception e) {
        return "Not Found / Error"
    }
}

// 実行前にツールバージョンを取得
def ver_fastp = get_version(params.fastp_path, "--version")
def ver_star  = get_version(params.star_path, "--version")
def ver_fc    = get_version(params.featurecounts_path, "-v")


// --- 入力チェック ---
if (!params.ref_gtf) {
    log.error "🚫 エラー: --ref_gtf（GTFファイル）は必須です！ncRNA検出のためにも必要です。"
    error "Missing GTF file."
}
if (!params.star_index && !params.ref_fasta) {
    log.error "🚫 エラー: インデックスを作成するための --ref_fasta が指定されていません。"
    error "Missing FASTA file."
}

// ログ出力
log.info """
          R N A - S E Q   F L E X I B L E   P I P E L I N E
          =================================================
          Pipeline Version : ${pipeline_version}
          =================================================
          SampleSheet   : ${params.samplesheet}
          FASTQ Dir     : ${params.fastq_dir}
          Ref GTF       : ${params.ref_gtf}
          Output Dir    : ${params.outdir}
          CPUs          : ${params.cpus}
          -------------------------------------------------
          [ Reference Info ]
          Ref Mode      : ${params.star_index ? "Use Existing Index" : "Build New Index from FASTA"}
          Ref FASTA     : ${params.ref_fasta ?: "N/A (Using Index)"}
          STAR Index    : ${params.star_index ?: "Create in: " + params.star_index_dir}
          
          [ Tools Versions ]
          fastp         : ${ver_fastp}
          STAR          : ${ver_star}
          featureCounts : ${ver_fc}
          =================================================
          """
          .stripIndent()

/*
 * ================================================================================================
 * ワークフロー定義
 * ================================================================================================
 */
workflow {
    // 1. サンプルシートの読み込み
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
        .set { ch_reads }

    // 2. FASTP実行
    FASTP( ch_reads )

    // 3. インデックスの準備（ロジック強化）
    def ch_index
    if (params.star_index) {
        // 既存インデックスを使う
        ch_index = Channel.fromPath(params.star_index).first()
    } else {
        // 新規作成（FASTAとGTFを使用）
        STAR_INDEX( params.ref_fasta, params.ref_gtf )
        ch_index = STAR_INDEX.out.index
    }

    // 4. STARアライメント
    // GTFをここでも渡すことで、ncRNA等のスプライシング精度を向上させます
    STAR_ALIGN( FASTP.out.reads, ch_index, params.ref_gtf )

    // 5. カウント処理
    // GTFに含まれるすべての遺伝子タイプ（protein_coding, lncRNA, miRNA等）をカウント
    FEATURECOUNTS( STAR_ALIGN.out.bam, params.ref_gtf )
}

/*
 * ================================================================================================
 * プロセス定義
 * ================================================================================================
 */

process FASTP {
    publishDir "${params.outdir}/fastp/${sample_id}", mode: 'copy'
    cpus params.cpus

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_Read{1,2}_trimmed.fastq.gz"), emit: reads
    path "${sample_id}.fastp.html", emit: html
    path "${sample_id}.fastp.json", emit: json

    script:
    def adapter_opt = params.adapter_fasta ? "--adapter_fasta ${params.adapter_fasta}" : ""
    
    """
    ${params.fastp_path} \\
        -i ${reads[0]} \\
        -I ${reads[1]} \\
        -o ${sample_id}_Read1_trimmed.fastq.gz \\
        -O ${sample_id}_Read2_trimmed.fastq.gz \\
        -h ${sample_id}.fastp.html \\
        -j ${sample_id}.fastp.json \\
        --thread ${task.cpus} \\
        --trim_poly_g \\
        --trim_poly_x \\
        ${adapter_opt}
    """
}

process STAR_INDEX {
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
    cpus params.cpus

    input:
    tuple val(sample_id), path(reads)
    path index
    path gtf  // ここでGTFを受け取るように追加！

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam
    path "${sample_id}_Log.final.out", emit: log

    script:
    // --sjdbGTFfile を追加して、マッピング時にGTF情報を活用するように変更
    """
    ${params.star_path} \\
        --genomeDir ${index} \\
        --readFilesIn ${reads[0]} ${reads[1]} \\
        --runThreadN ${task.cpus} \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${sample_id}_ \\
        --sjdbGTFfile ${gtf} \\
        --outSAMtype BAM SortedByCoordinate
    
    mv ${sample_id}_Aligned.sortedByCoord.out.bam ${sample_id}.bam
    """
}

process FEATURECOUNTS {
    publishDir "${params.outdir}/featurecounts", mode: 'copy'
    cpus params.cpus

    input:
    tuple val(sample_id), path(bam)
    path gtf

    output:
    path "${sample_id}_featurecounts.txt"
    path "${sample_id}_featurecounts.txt.summary"

    script:
    // -t exon -g gene_id はデフォルトですが、ncRNAもしっかり拾います
    """
    ${params.featurecounts_path} \\
        -p \\
        --countReadPairs \\
        -t exon \\
        -g gene_id \\
        -a ${gtf} \\
        -o ${sample_id}_featurecounts.txt \\
        -T ${task.cpus} \\
        ${bam}
    """
}