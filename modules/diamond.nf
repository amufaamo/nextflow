process DIAMOND {
    container 'quay.io/biocontainers/diamond:2.1.8--h43eeafb_0'
    publishDir "${params.outdir}/annotation/diamond", mode: 'copy'
    cpus params.cpus

    input:
    path pep_fasta
    path nr_db

    output:
    path "diamond_blastp.tsv", emit: tbl

    script:
    """
    diamond blastp -q ${pep_fasta} -d ${nr_db} -o diamond_blastp.tsv -f 6 -p ${task.cpus} --max-target-seqs 1 --evalue 1e-5
    """
}
