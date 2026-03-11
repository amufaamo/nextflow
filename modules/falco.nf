process FALCO_PRE {
    tag "${sample_id}"
    // FastQC/Falco用コンテナ
    container 'quay.io/biocontainers/falco:1.2.5--h077b44d_0'
    publishDir "${params.outdir}/falco_pre/${sample_id}", mode: 'copy'
    cpus params.cpus

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_report.html", emit: html
    path "*_data.txt", emit: txt

    script:
    def is_single = reads instanceof Path || reads.size() == 1
    def in_read1 = reads instanceof Path ? reads : reads[0]

    if (is_single) {
        """
        falco --outdir . -t ${task.cpus} ${in_read1}
        mv fastqc_report.html ${sample_id}_R1_pre_report.html
        mv fastqc_data.txt ${sample_id}_R1_pre_data.txt
        """
    } else {
        def in_read2 = reads[1]
        """
        falco --outdir . -t ${task.cpus} ${in_read1}
        mv fastqc_report.html ${sample_id}_R1_pre_report.html
        mv fastqc_data.txt ${sample_id}_R1_pre_data.txt

        falco --outdir . -t ${task.cpus} ${in_read2}
        mv fastqc_report.html ${sample_id}_R2_pre_report.html
        mv fastqc_data.txt ${sample_id}_R2_pre_data.txt
        """
    }
}

process FALCO_POST {
    tag "${sample_id}"
    // FastQC/Falco用コンテナ
    container 'quay.io/biocontainers/falco:1.2.5--h077b44d_0'
    publishDir "${params.outdir}/falco_post/${sample_id}", mode: 'copy'
    cpus params.cpus

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_report.html", emit: html
    path "*_data.txt", emit: txt

    script:
    def is_single = reads instanceof Path || reads.size() == 1
    def in_read1 = reads instanceof Path ? reads : reads[0]

    if (is_single) {
        """
        falco --outdir . -t ${task.cpus} ${in_read1}
        mv fastqc_report.html ${sample_id}_R1_post_report.html
        mv fastqc_data.txt ${sample_id}_R1_post_data.txt
        """
    } else {
        def in_read2 = reads[1]
        """
        falco --outdir . -t ${task.cpus} ${in_read1}
        mv fastqc_report.html ${sample_id}_R1_post_report.html
        mv fastqc_data.txt ${sample_id}_R1_post_data.txt

        falco --outdir . -t ${task.cpus} ${in_read2}
        mv fastqc_report.html ${sample_id}_R2_post_report.html
        mv fastqc_data.txt ${sample_id}_R2_post_data.txt
        """
    }
}
