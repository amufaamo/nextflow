#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * ================================================================================================
 * 🌟 Universal RNA-Seq Pipeline v2.0 (Ref-based & De novo) 🌟
 * ================================================================================================
 */

def pipeline_version = "2.0.0"

// --- ツールバージョン定義 ---
def fastp_version   = "0.23.4"
def star_version    = "2.7.10b"
def subread_version = "2.0.1"
def trinity_version = "2.15.1"
def busco_version   = "5.5.0"
def salmon_version  = "1.10.1"

// --- Apptainer チェック ---
if (workflow.profile == 'standard' || workflow.profile == 'apptainer' || workflow.profile == 'singularity') {
    try {
        def proc = "apptainer --version".execute()
        proc.waitFor()
        if (proc.exitValue() != 0) throw new Exception()
    } catch (Exception e) {
        log.error "🚫 エラー: Apptainerが見つかりません！"
        System.exit(1)
    }
}

// --- モジュールのインポート ---
include { FALCO_PRE; FALCO_POST } from '../modules/falco.nf'
include { FASTP } from '../modules/fastp.nf'
include { SUMMARY } from '../modules/summary.nf'
include { STAR_INDEX; STAR_ALIGN } from '../modules/star.nf'
include { FEATURECOUNTS } from '../modules/featurecounts.nf'
include { TRINITY } from '../modules/trinity.nf'
include { BUSCO } from '../modules/busco.nf'
include { SALMON_INDEX; SALMON_QUANT } from '../modules/salmon.nf'
include { CORSET } from '../modules/corset.nf'
include { DESEQ2 } from '../modules/deseq2.nf'
include { TRANSDECODER } from '../modules/transdecoder.nf'
include { EGGNOG_MAPPER } from '../modules/eggnog.nf'
include { DIAMOND } from '../modules/diamond.nf'
include { INTERPROSCAN } from '../modules/interproscan.nf'

// --- パラメータ設定 (configで上書きされます) ---
params.samplesheet     = './samples.csv'
params.outdir          = "results"
params.denovo          = false 
params.read_type       = 'paired' // 'paired', 'read1', 'read2'

// --- De novo Annotation 用 データベース群 ---
params.eggnog_db       = null
params.nr_db           = null
params.interpro_db     = null

// --- 入力チェック ---
if (params.denovo) {
    log.info "🚀 Mode: De novo Assembly (Trinity -> BUSCO -> Salmon)"
} else {
    log.info "🚀 Mode: Reference-based (STAR -> featureCounts)"
}

log.info """
          R N A - S E Q   P I P E L I N E  (v${pipeline_version})
          =================================================
          [Mode]
          De novo Assembly : ${params.denovo}

          [Run Info]
          SampleSheet      : ${params.samplesheet}
          Output Dir       : ${params.outdir}
          =================================================
          """
          .stripIndent()

/*
 * ================================================================================================
 * ワークフロー定義
 * ================================================================================================
 */
workflow {
    // 1. CSV読み込み
    Channel
        .fromPath( params.samplesheet )
        .splitCsv( header:true )
        .map { row ->
            def sample_id = row.sample
            def reads
            def read_type_opt = params.read_type ?: (params.single_end ? 'read1' : 'paired')
            if (read_type_opt == 'read1') {
                reads = [ file("${params.fastq_dir}/${row.fastq_1}") ]
            } else if (read_type_opt == 'read2') {
                reads = [ file("${params.fastq_dir}/${row.fastq_2}") ]
            } else {
                reads = [ file("${params.fastq_dir}/${row.fastq_1}"), file("${params.fastq_dir}/${row.fastq_2}") ]
            }
            tuple( sample_id, reads )
        }
        .set { ch_reads }

    // 2. フィルタリング前の QC (Falco)
    FALCO_PRE( ch_reads )

    // 3. QC & Trimming (fastp)
    FASTP( ch_reads )

    // 4. フィルタリング後の QC (Falco)
    FALCO_POST( FASTP.out.reads )

    if (params.denovo) {
        // ==========================================
        // 🌿 De novo Route
        // ==========================================
        
        // ★修正ポイント: .set{} をやめて、変数へ直接代入(=)に変更
        def ch_r1_list
        def ch_r2_list

        def is_single = params.single_end || params.read_type == 'read1' || params.read_type == 'read2'

        if (is_single) {
            // Single Endの場合
            ch_r1_list = FASTP.out.reads
                .map { id, files -> files instanceof List ? files[0] : files }
                .collect()
            
            // R2は空のリストを渡す
            ch_r2_list = Channel.value( [] )
        } else {
            // Paired Endの場合
            // R1ファイルを抽出してリスト化
            ch_r1_list = FASTP.out.reads
                .map { id, files -> files[0] } // R1を取り出す
                .collect() 

            // R2ファイルを抽出してリスト化
            ch_r2_list = FASTP.out.reads
                .map { id, files -> files[1] } // R2を取り出す
                .collect() 
        }

        // Trinity実行 (安全なリストを渡す)
        TRINITY( ch_r1_list, ch_r2_list )
        
        // アセンブリの品質評価
        BUSCO( TRINITY.out.fasta, params.busco_lineage )
        
        // ==========================================
        // 🔀 【分岐スタート】Nextflowが自動で同時並行で走らせます
        // ==========================================
        
        // 🔴 【定量・クラスタリングルート】
        // 定量ルートでは、TranscriptレベルのFASTA（ Trinityの出力 ）を使います
        SALMON_INDEX( TRINITY.out.fasta )
        SALMON_QUANT( FASTP.out.reads, SALMON_INDEX.out.index )
        
        // Corsetで遺伝子レベルにクラスタリング
        CORSET( SALMON_QUANT.out.eq_classes.collect() )
        
        // DESeq2による発現変動解析
        DESEQ2( CORSET.out.counts )

        // 🔵 【アノテーションルート】
        // FASTAからペプチド(ORF)配列を抽出
        TRANSDECODER( TRINITY.out.fasta )
        
        // 以下、外部DBのパスが指定されている場合のみ各アノテーションツールを実行
        if (params.eggnog_db) {
            EGGNOG_MAPPER( TRANSDECODER.out.pep, params.eggnog_db )
        }
        
        if (params.nr_db) {
            DIAMOND( TRANSDECODER.out.pep, params.nr_db )
        }
        
        if (params.interpro_db) {
            INTERPROSCAN( TRANSDECODER.out.pep, params.interpro_db )
        }

    } else {
        // ==========================================
        // 🗺️ Reference-based Route
        // ==========================================
        def ch_index
        if (params.star_index) {
            ch_index = Channel.fromPath(params.star_index).first()
        } else {
            STAR_INDEX( params.ref_fasta, params.ref_gtf )
            ch_index = STAR_INDEX.out.index
        }
        STAR_ALIGN( FASTP.out.reads, ch_index, params.ref_gtf )
        FEATURECOUNTS( STAR_ALIGN.out.bam, params.ref_gtf )

        // Summary出力
        ch_falco_pre  = FALCO_PRE.out.txt.collect()
        ch_falco_post = FALCO_POST.out.txt.collect()
        ch_fastp      = FASTP.out.json.collect()
        ch_star       = STAR_ALIGN.out.log.collect()
        ch_fc         = FEATURECOUNTS.out.summary.collect()
        
        SUMMARY( ch_falco_pre, ch_falco_post, ch_fastp, ch_star, ch_fc )
    }
}