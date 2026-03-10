process TRANSDECODER {
    container 'quay.io/biocontainers/transdecoder:5.7.1--pl5321hdfd78af_0'
    publishDir "${params.outdir}/annotation/transdecoder", mode: 'copy'
    cpus params.cpus

    input:
    path transcripts_fasta

    output:
    path "*.pep", emit: pep
    path "*.cds", emit: cds

    script:
    """
    TransDecoder.LongOrfs -t ${transcripts_fasta}
    TransDecoder.Predict -t ${transcripts_fasta} --single_best_only
    """
}
