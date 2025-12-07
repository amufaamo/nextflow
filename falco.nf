#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//params.fastq = "/home/masakazu/fastq*.fastq.gz"
params.fastq = "*.f*q.gz"
params.fastq_dir = "/home/masakazu/fastq"
params.outdir = "/home/masakazu/falco/output"
params.falco_path = "/home/masakazu/miniconda3/envs/falco/bin/falco"

process runFalco {
    tag { fastq.baseName }
    publishDir "${params.outdir}/", mode: 'copy'

    input:
    path fastq

    output:
    // Falcoが実際に出力するファイル名に修正
    path "${fastq.baseName.replaceAll('\\.fastq$|\\.fq$', '')}_report.html"
    path "${fastq.baseName.replaceAll('\\.fastq$|\\.fq$', '')}_data.txt"

    script:
    def base = fastq.baseName.replaceAll('\\.fastq$|\\.fq$', '')
    """
    ${params.falco_path} --outdir . ${fastq}
    mv fastqc_report.html ${base}_report.html
    mv fastqc_data.txt ${base}_data.txt
    """
}

workflow {
    def input_pattern = "${params.fastq_dir}/${params.fastq}"
    log.info "Using input file pattern: ${input_pattern}"
    fastq_ch = Channel.fromPath(input_pattern, checkIfExists: true)
    runFalco(fastq_ch)
}