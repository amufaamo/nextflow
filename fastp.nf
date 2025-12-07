#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// --- パラメータ定義 (デフォルト値) ---
// コマンドラインで --input_dir path/to/fastq のように上書き可
params.input_dir = System.getenv('HOME') + '/fastq'      // 入力fastqファイルが含まれるディレクトリ
params.output_dir = System.getenv('HOME') + '/fastq'   // 出力ディレクトリ
params.sample_id = false
params.read1 = false
params.read2 = false
params.adapter_fasta = System.getenv('HOME') + '/fasta/adapter.fasta' // アダプターFASTAファイル
params.fastp_path = System.getenv('HOME') + '/miniconda3/envs/fastp/bin/fastp' // fastp実行ファイルのパス

// --- パラメータのチェック ---
// --read1 と --read2 が指定されているか確認
if (!params.read1 || !params.read2) {
    // error() でパイプラインを停止し、メッセージを表示
    error "エラー: 入力ファイルを --read1 と --read2 で指定してください。\n例: nextflow run main.nf --read1 R1.fastq.gz --read2 R2.fastq.gz"
}

// --- 入力ディレクトリのパスオブジェクトを作成し、存在確認 ---
def input_dir_path = file(params.input_dir) // file() で Path オブジェクトに変換
if (!input_dir_path.isDirectory()) { // ディレクトリかどうかチェック
    error "エラー: 指定された入力ディレクトリが見つからないか、ディレクトリではありません: ${params.input_dir}"
}

// --- 完全なファイルパスを構築 ---
// input_dir_path (ディレクトリのPathオブジェクト) とファイル名を結合
def read1_full_path = input_dir_path.resolve(params.read1)
def read2_full_path = input_dir_path.resolve(params.read2)

// --- ファイルオブジェクトを作成 ---
// 完全なパスからファイルオブジェクトを作成
def read1_file = file(read1_full_path)
def read2_file = file(read2_full_path)

// --- ファイルの存在チェック ---
// 完全なパスで存在を確認
if (!read1_file.exists()) {
    error "エラー: Read1ファイルが見つかりません: ${read1_full_path}"
}
if (!read2_file.exists()) {
    error "エラー: Read2ファイルが見つかりません: ${read2_full_path}"
}


// --- サンプルIDを決定 ---
def final_sample_id
if (params.sample_id) {
    final_sample_id = params.sample_id
} else {
    // Read1のファイル名 (params.read1) から拡張子やペア識別子を除去
    final_sample_id = params.read1.replaceAll(/.fastq.gz$|.fq.gz$|.fastq$|.fq$/, '').replaceAll(/_R1$|_1$|.1$/, '')
    if (final_sample_id == params.read1.replaceAll(/.gz$/, '')) { // うまく除去できなかった場合
         log.warn "サンプルIDをファイル名 '${params.read1}' から自動特定できませんでした。ファイル名（拡張子除去）を使用します: ${final_sample_id}"
    }
}




println """
         F A S T P   P I P E L I N E (Simple)
         ====================================
         Input Directory : ${params.input_dir}
         Output Directory: ${params.output_dir}
         Adapter FASTA   : ${params.adapter_fasta}
         FASTP Path      : ${params.fastp_path}
         """



// --- 入力チャンネルの作成 (単一要素) ---
// read1_file, read2_file は完全なパスを指すファイルオブジェクト
def reads_ch = Channel.of([ final_sample_id, [read1_file, read2_file] ])


// --- FASTP プロセス定義 ---
process runFastp {
    tag "$sample_id"
    publishDir "${params.output_dir}", mode: 'copy' // 出力モードを直接指定

    input:
    tuple val(sample_id), path(reads)
    path adapter_fasta

    output:
    tuple val(sample_id), path("${sample_id}_R1_trimmed.fastq.gz"), path("${sample_id}_R2_trimmed.fastq.gz"), emit: trimmed_reads
    path "${sample_id}.fastp.json", emit: json_report
    path "${sample_id}.fastp.html", emit: html_report

    script:
    def (r1, r2) = reads
    def output_r1 = "${sample_id}_R1_trimmed.fastq.gz"
    def output_r2 = "${sample_id}_R2_trimmed.fastq.gz"
    def json_out = "${sample_id}.fastp.json"
    def html_out = "${sample_id}.fastp.html"

    """
    echo "Processing sample: ${sample_id}"
    ${params.fastp_path} \\
        -i ${r1} \\
        -I ${r2} \\
        -o ${output_r1} \\
        -O ${output_r2} \\
        --trim_poly_g \\
        --trim_poly_x \\
        --adapter_fasta ${adapter_fasta} \\
        -j ${json_out} \\
        -h ${html_out}
    echo "Finished sample: ${sample_id}"
    """
}

// --- アダプターファイルのチャンネルを作成 --- ★★★この行があるか確認★★★
adapter_fasta_ch = file(params.adapter_fasta) // この行が workflow の前に必要
// --- ワークフロー定義 ---


 workflow {
     println "Debug: reads_ch content:" // デバッグが不要ならこの行も削除してOK
     reads_ch.view()                  // デバッグが不要ならこの行も削除してOK

     // FASTP プロセス呼び出し (adapter_fasta_ch はそのまま渡す)
    runFastp(reads_ch, adapter_fasta_ch) // ★★★ プロセス名を "runFastp" に修正 ★★★
 }