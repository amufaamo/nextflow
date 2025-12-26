#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * ================================================================================================
 * 🧬 Single Cell RNA-Seq Pipeline (Pipseeker + Seurat)
 * ================================================================================================
 */

// モジュールの読み込み (../module/ を参照)
include { PIPSEEKER_INDEX; PIPSEEKER_FULL } from '../module/pipseeker.nf'
include { SEURAT_PROCESS }                  from '../module/seurat.nf'

// 必須パラメータのチェック
if (!params.fasta || !params.gtf) {
    log.error "❌ エラー: --fasta と --gtf は必須です！"
    System.exit(1)
}

log.info """
          S C - R N A - S E Q   P I P E L I N E
          =================================================
          Ref Fasta    : ${params.fasta}
          Ref GTF      : ${params.gtf}
          Outdir       : ${params.outdir}
          Chemistry    : ${params.sc_chemistry}
          =================================================
          """.stripIndent()

workflow {
    
    // =================================================================
    // 1. 入力データの読み込み (モード分岐)
    // =================================================================
    def ch_reads

    if (params.r1 && params.r2) {
        // 【A】 単一サンプルモード (コマンド引数でファイルを直接指定)
        log.info "�� Mode: Single Sample Input (--r1, --r2)"
        
        ch_reads = Channel.of(
            tuple( params.sample_name, [ file(params.r1), file(params.r2) ] )
        )

    } else if (params.samplesheet && params.fastq_dir) {
        // 【B】 サンプルシートモード (CSV + ディレクトリ指定)
        log.info "🚀 Mode: Batch Processing (--samplesheet)"
        
        ch_reads = Channel
            .fromPath( params.samplesheet )
            .splitCsv( header:true )
            .map { row ->
                tuple( row.sample, [ file("${params.fastq_dir}/${row.fastq_1}"), file("${params.fastq_dir}/${row.fastq_2}") ] )
            }

    } else {
        log.error "❌ エラー: 入力データが指定されていません。"
        log.error "   単一サンプル: --r1 <file> --r2 <file>"
        log.error "   複数サンプル: --samplesheet <csv> --fastq_dir <dir>"
        System.exit(1)
    }

    // =================================================================
    // 2. リファレンスの準備
    // =================================================================
    // インデックス作成 (Pipseekerはbuildmaprefが必要)
    // ※今回は毎回指定されたFASTA/GTFから作成するフローにしています
    PIPSEEKER_INDEX( params.fasta, params.gtf )
    
    // =================================================================
    // 3. 解析実行
    // =================================================================
    // Pipseeker
    PIPSEEKER_FULL( ch_reads, PIPSEEKER_INDEX.out.index )

    // Seurat
    SEURAT_PROCESS( PIPSEEKER_FULL.out.out_dir )
}
