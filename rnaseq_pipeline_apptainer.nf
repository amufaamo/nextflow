#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * ================================================================================================
 * 🌟 Universal RNA-Seq Pipeline v1.4.1 (Fix: Tool Paths) 🌟
 * ================================================================================================
 */

def pipeline_version = "1.4.2"

// --- 必須パラメータ ---
params.ref_gtf     = null

// --- インデックス関連 ---
params.ref_fasta   = null
params.star_index  = null

// --- オプションパラメータ ---
params.samplesheet = './samples.csv'
params.fastq_dir   = "${System.getProperty("user.home")}/fastq"
params.outdir      = "results"
params.cpus        = 8
params.single_end  = false 

// featureCountsの集計単位 (デフォルト: gene_id)
params.fc_group_features = 'gene_id'

// --- ツールコマンド名 (★ここが重要！これが無いとエラーになります) ---
params.fastp_path         = 'fastp'
params.star_path          = 'STAR'
params.featurecounts_path = 'featureCounts'

params.adapter_fasta      = "${System.getProperty("user.home")}/fasta/adapter.fasta"
params.star_index_dir     = "./star_index_new" 

// --- 入力チェック ---
if (!params.ref_gtf) {
    log.error "🚫 エラー: --ref_gtf（GTFファイル）は必須です！"
    error "Missing GTF file."
}
if (!params.star_index && !params.ref_fasta) {
    log.error "🚫 エラー: --ref_fasta または --star_index が必要です。"
    error "Missing FASTA file."
}

// ★実行時にロゴとバージョンを表示
log.info """
          R N A - S E Q   P I P E L I N E 
          =================================================
          Pipeline Version : ${pipeline_version}
          SampleSheet      : ${params.samplesheet}
          Output Dir       : ${params.outdir}
          Single End       : ${params.single_end}
          Feature ID (-g)  : ${params.fc_group_features}
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
    tuple val(sample_id), path("${sample_id}*.fastq.gz"), emit: reads
    path "${sample_id}.fastp.html", emit: html
    path "${sample_id}.fastp.json", emit: json

    script:
    def adapter_opt = params.adapter_fasta ? "--adapter_fasta ${params.adapter_fasta}" : ""
    
    if (params.single_end) {
        """
        ${params.fastp_path} -i ${reads[0]} -o ${sample_id}_trimmed.fastq.gz -h ${sample_id}.fastp.html -j ${sample_id}.fastp.json --thread ${task.cpus} --trim_poly_g --trim_poly_x ${adapter_opt}
        """
    } else {
        """
        ${params.fastp_path} -i ${reads[0]} -I ${reads[1]} -o ${sample_id}_Read1_trimmed.fastq.gz -O ${sample_id}_Read2_trimmed.fastq.gz -h ${sample_id}.fastp.html -j ${sample_id}.fastp.json --thread ${task.cpus} --trim_poly_g --trim_poly_x ${adapter_opt}
        """
    }
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
    ${params.star_path} --runMode genomeGenerate --genomeDir index --genomeFastaFiles ${fasta} --sjdbGTFfile ${gtf} --runThreadN ${task.cpus}
    """
}

process STAR_ALIGN {
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
    """
    ${params.star_path} --genomeDir ${index} --readFilesIn ${input_reads} --runThreadN ${task.cpus} --readFilesCommand zcat --outFileNamePrefix ${sample_id}_ --sjdbGTFfile ${gtf} --outSAMtype BAM SortedByCoordinate
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
    def paired_opts = params.single_end ? "" : "-p --countReadPairs"

    """
    ${params.featurecounts_path} \
        ${paired_opts} \
        -t exon \
        -g ${params.fc_group_features} \
        -a ${gtf} \
        -o ${sample_id}_featurecounts.txt \
        -T ${task.cpus} \
        ${bam}
    """
}