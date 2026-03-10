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
params.falco_path         = '/home/masakazu/miniconda3/envs/falco/bin/falco'
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
def ver_falco = get_version(params.falco_path, "-v")


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
          falco         : ${ver_falco}
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

    // 2. フィルタリング前のQC (Falco)
    FALCO_PRE( ch_reads )

    // 3. FASTP実行
    FASTP( ch_reads )

    // 4. フィルタリング後のQC (Falco)
    FALCO_POST( FASTP.out.reads )

    // 5. インデックスの準備（ロジック強化）
    def ch_index
    if (params.star_index) {
        // 既存インデックスを使う
        ch_index = Channel.fromPath(params.star_index).first()
    } else {
        // 新規作成（FASTAとGTFを使用）
        STAR_INDEX( params.ref_fasta, params.ref_gtf )
        ch_index = STAR_INDEX.out.index
    }

    // 6. STARアライメント
    // GTFをここでも渡すことで、ncRNA等のスプライシング精度を向上させます
    STAR_ALIGN( FASTP.out.reads, ch_index, params.ref_gtf )

    // 7. カウント処理
    // GTFに含まれるすべての遺伝子タイプ（protein_coding, lncRNA, miRNA等）をカウント
    FEATURECOUNTS( STAR_ALIGN.out.bam, params.ref_gtf )

    // 8. ログのサマリー作成（TSVに集計）
    ch_falco_pre  = FALCO_PRE.out.txt.collect()
    ch_falco_post = FALCO_POST.out.txt.collect()
    ch_fastp      = FASTP.out.json.collect()
    ch_star       = STAR_ALIGN.out.log.collect()
    ch_fc         = FEATURECOUNTS.out.summary.collect()

    SUMMARY( ch_falco_pre, ch_falco_post, ch_fastp, ch_star, ch_fc )
}

/*
 * ================================================================================================
 * プロセス定義
 * ================================================================================================
 */

process FALCO_PRE {
    tag "${sample_id}"
    publishDir "${params.outdir}/falco_pre/${sample_id}", mode: 'copy'
    cpus params.cpus

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_report.html", emit: html
    path "*_data.txt", emit: txt

    script:
    """
    ${params.falco_path} --outdir . -t ${task.cpus} ${reads[0]}
    mv fastqc_report.html ${sample_id}_R1_pre_report.html
    mv fastqc_data.txt ${sample_id}_R1_pre_data.txt

    ${params.falco_path} --outdir . -t ${task.cpus} ${reads[1]}
    mv fastqc_report.html ${sample_id}_R2_pre_report.html
    mv fastqc_data.txt ${sample_id}_R2_pre_data.txt
    """
}

process FALCO_POST {
    tag "${sample_id}"
    publishDir "${params.outdir}/falco_post/${sample_id}", mode: 'copy'
    cpus params.cpus

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_report.html", emit: html
    path "*_data.txt", emit: txt

    script:
    """
    ${params.falco_path} --outdir . -t ${task.cpus} ${reads[0]}
    mv fastqc_report.html ${sample_id}_R1_post_report.html
    mv fastqc_data.txt ${sample_id}_R1_post_data.txt

    ${params.falco_path} --outdir . -t ${task.cpus} ${reads[1]}
    mv fastqc_report.html ${sample_id}_R2_post_report.html
    mv fastqc_data.txt ${sample_id}_R2_post_data.txt
    """
}

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
    path "${sample_id}_featurecounts.txt.summary", emit: summary

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

process SUMMARY {
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    path falco_pre_txt
    path falco_post_txt
    path fastp_json
    path star_logs
    path fc_summaries

    output:
    path "summary.tsv"

    script:
    """
    #!/usr/bin/env python3
    import glob
    
    data = {}
    samples = set()
    TAB = chr(9)
    NL = chr(10)
    
    def add_metric(sample, metric, value):
        if sample not in data:
            data[sample] = {}
        data[sample][metric] = value
        samples.add(sample)
    
    for f in glob.glob("*_R1_pre_data.txt"):
        sample = f.replace("_R1_pre_data.txt", "")
        with open(f) as fh:
            for line in fh:
                if line.startswith("Total Sequences"):
                    add_metric(sample, "1_Falco_Pre_Total_Seq", line.strip().split(TAB)[1])
                elif line.startswith("%GC"):
                    add_metric(sample, "1_Falco_Pre_GC_Percent", line.strip().split(TAB)[1])
    
    for f in glob.glob("*_R1_post_data.txt"):
        sample = f.replace("_R1_post_data.txt", "")
        with open(f) as fh:
            for line in fh:
                if line.startswith("Total Sequences"):
                    add_metric(sample, "2_Falco_Post_Total_Seq", line.strip().split(TAB)[1])
                elif line.startswith("%GC"):
                    add_metric(sample, "2_Falco_Post_GC_Percent", line.strip().split(TAB)[1])
    
    for f in glob.glob("*_Log.final.out"):
        sample = f.replace("_Log.final.out", "")
        with open(f) as fh:
            for line in fh:
                if "Number of input reads" in line:
                    add_metric(sample, "3_STAR_Input_Reads", line.split("|")[1].strip())
                elif "Uniquely mapped reads number" in line:
                    add_metric(sample, "3_STAR_Uniquely_Mapped", line.split("|")[1].strip())
                elif "Uniquely mapped reads %" in line:
                    add_metric(sample, "3_STAR_Uniquely_Mapped_Percent", line.split("|")[1].strip())
    
    for f in glob.glob("*_featurecounts.txt.summary"):
        sample = f.replace("_featurecounts.txt.summary", "")
        with open(f) as fh:
            for i, line in enumerate(fh):
                if i == 0: continue
                parts = line.strip().split(TAB)
                if len(parts) >= 2:
                    add_metric(sample, "4_FC_" + parts[0], parts[1])
    
    with open("summary.tsv", "w") as out:
        samples = sorted(list(samples))
        out.write("Metric" + TAB + TAB.join(samples) + NL)
        seen_metrics = set()
        for s in data:
            seen_metrics.update(data[s].keys())
        for m in sorted(list(seen_metrics)):
            row = [m] + [str(data.get(s, {}).get(m, "NA")) for s in samples]
            out.write(TAB.join(row) + NL)
    """
}