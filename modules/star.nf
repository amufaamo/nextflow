process STAR_INDEX {
    // STAR container definition
    container 'quay.io/biocontainers/star:2.7.10b--h6b7c446_1'

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
    STAR --runMode genomeGenerate \
        --genomeDir index \
        --genomeFastaFiles ${fasta} \
        --sjdbGTFfile ${gtf} \
        --runThreadN ${task.cpus} \
        --genomeSAindexNbases ${params.star_index_nbases}
    """
}

process STAR_ALIGN {
    tag "$sample_id"
    // STAR container definition
    container 'quay.io/biocontainers/star:2.7.10b--h6b7c446_1'

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
    STAR --genomeDir ${index} \
        --readFilesIn ${input_reads} \
        --runThreadN ${task.cpus} \
        --readFilesCommand zcat \
        --outFileNamePrefix ${sample_id}_ \
        --sjdbGTFfile ${gtf} \
        --outSAMtype BAM SortedByCoordinate
    
    mv ${sample_id}_Aligned.sortedByCoord.out.bam ${sample_id}.bam
    """
}
