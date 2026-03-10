process INTERPROSCAN {
    container 'interpro/interproscan:5.63-95.0'
    publishDir "${params.outdir}/annotation/interproscan", mode: 'copy'
    cpus params.cpus

    input:
    path pep_fasta
    path interpro_db_dir

    output:
    path "*.tsv", emit: tsv
    path "*.gff3", emit: gff3

    script:
    """
    # Asterisks sometimes break InterProScan, so we sanitize the FASTA
    sed 's/\\*//g' ${pep_fasta} > clean.pep
    
    interproscan.sh -i clean.pep -f TSV,GFF3 -o interproscan_out -dp -cpu ${task.cpus} -d ${interpro_db_dir}
    """
}
