#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * ================================================================================================
 * 🌟 Universal RNA-Seq Pipeline v1.5.0
 * ================================================================================================
 */

def pipeline_version = "1.5.0"

// --- 使用するツールのバージョン定義（ログ表示用） ---
// ※実体は各モジュール内のコンテナタグと一致させてください
def fastp_version   = "0.23.4"
def star_version    = "2.7.10b"
def subread_version = "2.0.1"

// --- Apptainer チェック ---
// Apptainer/Singularityプロファイルが有効な場合、コマンドの有無を確認
if (workflow.profile == 'standard' || workflow.profile == 'apptainer' || workflow.profile == 'singularity') {
    try {
        def proc = "apptainer --version".execute()
        proc.waitFor()
        if (proc.exitValue() != 0) throw new Exception()
    } catch (Exception e) {
        log.error """
        ================================================================
        �� エラー: Apptainerが見つかりません！
        
        'apptainer' コマンドがインストールされているか、パスが通っているか確認してください。
        このパイプラインはデフォルトでApptainerを使用します。
        ================================================================
        """.stripIndent()
        System.exit(1)
    }
}

// --- モジュールのインポート（相対パス：一つ上の階層のmodulesを参照） ---
include { FASTP } from '../modules/fastp.nf'
include { STAR_INDEX; STAR_ALIGN } from '../modules/star.nf'
include { FEATURECOUNTS } from '../modules/featurecounts.nf'

// --- パラメータ設定 ---
params.ref_gtf         = null
params.ref_fasta       = null
params.star_index      = null
params.samplesheet     = './samples.csv'
params.fastq_dir       = "${System.getProperty("user.home")}/fastq"
params.outdir          = "results"
params.cpus            = 8
params.single_end      = false 
params.fc_group_features = 'gene_id'
params.adapter_fasta   = "${System.getProperty("user.home")}/fasta/adapter.fasta"
params.star_index_dir  = "./star_index_new" 

// --- 入力チェック ---
if (!params.ref_gtf) {
    log.error "🚫 エラー: --ref_gtf（GTFファイル）は必須です！"
    error "Missing GTF file."
}
if (!params.star_index && !params.ref_fasta) {
    log.error "🚫 エラー: --ref_fasta または --star_index が必要です。"
    error "Missing FASTA file."
}

// --- 実行時のロゴと情報表示 ---
log.info """
          R N A - S E Q   P I P E L I N E 
          =================================================
          Pipeline Version : ${pipeline_version}
          
          [Tools Version]
          fastp            : ${fastp_version}
          STAR             : ${star_version}
          featureCounts    : ${subread_version} (Subread)

          [Run Info]
          SampleSheet      : ${params.samplesheet}
          Output Dir       : ${params.outdir}
          Single End       : ${params.single_end}
          Feature ID (-g)  : ${params.fc_group_features}
          Container Engine : ${workflow.containerEngine ?: 'local'}
          =================================================
          """
          .stripIndent()

/*
 * ================================================================================================
 * ワークフロー定義
 * ================================================================================================
 */
workflow {
    Channel
        .fromPath( params.samplesheet )
        .splitCsv( header:true )
        .map { row ->
            def sample_id = row.sample
            def reads
            if (params.single_end) {
                reads = [ file("${params.fastq_dir}/${row.fastq_1}") ]
            } else {
                reads = [ file("${params.fastq_dir}/${row.fastq_1}"), file("${params.fastq_dir}/${row.fastq_2}") ]
            }
            tuple( sample_id, reads )
        }
        .set { ch_reads }

    FASTP( ch_reads )

    def ch_index
    if (params.star_index) {
        ch_index = Channel.fromPath(params.star_index).first()
    } else {
        STAR_INDEX( params.ref_fasta, params.ref_gtf )
        ch_index = STAR_INDEX.out.index
    }

    STAR_ALIGN( FASTP.out.reads, ch_index, params.ref_gtf )
    FEATURECOUNTS( STAR_ALIGN.out.bam, params.ref_gtf )
}
