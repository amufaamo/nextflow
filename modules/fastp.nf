process FASTP {
    tag "$sample_id"
    // Fastp container definition
    container 'quay.io/biocontainers/fastp:1.0.1--heae3180_0'
    
    
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
    
    def single_end = reads instanceof Path || reads.size() == 1
    def in_read1 = reads instanceof Path ? reads : reads[0]

    if (single_end) {
        """
        fastp -i ${in_read1} -o ${sample_id}_trimmed.fastq.gz \\
            -h ${sample_id}.fastp.html -j ${sample_id}.fastp.json \\
            --thread ${task.cpus} \\
            --trim_poly_g --trim_poly_x \\
            ${adapter_opt}
        """
    } else {
        def in_read2 = reads[1]
        """
        fastp -i ${in_read1} -I ${in_read2} \\
            -o ${sample_id}_Read1_trimmed.fastq.gz -O ${sample_id}_Read2_trimmed.fastq.gz \\
            -h ${sample_id}.fastp.html -j ${sample_id}.fastp.json \\
            --thread ${task.cpus} \\
            --trim_poly_g --trim_poly_x \\
            ${adapter_opt}
        """
    }
}
