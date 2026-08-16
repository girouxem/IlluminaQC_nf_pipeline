// modules/multiqc.nf
// Contains: MULTIQC

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    cpus   { params.global_cpus }
    memory { params.global_memory }
    
    input:
    path all_reports
    path multiqc_config
    
    output:
    path "multiqc_report.html"
    
    script:
    """
    multiqc . -o ./ -c ${multiqc_config}
    """
}
