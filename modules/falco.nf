process FALCO_PRE {
    tag "${sample_id}"
    // FastQC/Falco用コンテナ
    container 'quay.io/biocontainers/falco:1.2.1--hec16e2b_1'
    publishDir "${params.outdir}/falco_pre/${sample_id}", mode: 'copy'
    cpus params.cpus

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_report.html", emit: html
    path "*_data.txt", emit: txt

    script:
    if (params.single_end) {
        """
        falco --outdir . -t ${task.cpus} ${reads[0]}
        mv fastqc_report.html ${sample_id}_R1_pre_report.html
        mv fastqc_data.txt ${sample_id}_R1_pre_data.txt
        """
    } else {
        """
        falco --outdir . -t ${task.cpus} ${reads[0]}
        mv fastqc_report.html ${sample_id}_R1_pre_report.html
        mv fastqc_data.txt ${sample_id}_R1_pre_data.txt

        falco --outdir . -t ${task.cpus} ${reads[1]}
        mv fastqc_report.html ${sample_id}_R2_pre_report.html
        mv fastqc_data.txt ${sample_id}_R2_pre_data.txt
        """
    }
}

process FALCO_POST {
    tag "${sample_id}"
    // FastQC/Falco用コンテナ
    container 'quay.io/biocontainers/falco:1.2.1--hec16e2b_1'
    publishDir "${params.outdir}/falco_post/${sample_id}", mode: 'copy'
    cpus params.cpus

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_report.html", emit: html
    path "*_data.txt", emit: txt

    script:
    if (params.single_end) {
        """
        falco --outdir . -t ${task.cpus} ${reads[0]}
        mv fastqc_report.html ${sample_id}_R1_post_report.html
        mv fastqc_data.txt ${sample_id}_R1_post_data.txt
        """
    } else {
        """
        falco --outdir . -t ${task.cpus} ${reads[0]}
        mv fastqc_report.html ${sample_id}_R1_post_report.html
        mv fastqc_data.txt ${sample_id}_R1_post_data.txt

        falco --outdir . -t ${task.cpus} ${reads[1]}
        mv fastqc_report.html ${sample_id}_R2_post_report.html
        mv fastqc_data.txt ${sample_id}_R2_post_data.txt
        """
    }
}
