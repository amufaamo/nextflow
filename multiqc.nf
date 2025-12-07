nextflow.enable.dsl=2

params.output_dir = 'falco/output'
params.multiqc_report = 'multiqc_report'

process RUN_MULTIQC {
    tag "Running MultiQC on ${params.output_dir}"

    input:
    path input_dir from params.output_dir

    output:
    path "${params.multiqc_report}.html"
    path "${params.multiqc_report}_data"

    script:
    """
    /home/masakazu/miniconda3/envs/multiqc/bin/multiqc \
        --force \
        --outdir . \
        --filename ${params.multiqc_report} \
        ${input_dir}
    """
}

workflow {
    RUN_MULTIQC()
}
