#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Case‑sensitive match first (we add case-insensitive after baseline works)
params.reads  = "data/*_{r1,r2}.fastq.gz"
params.outdir = "${baseDir}/results"

process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    input:
    tuple val(sample_id), path(reads)
    output:
    path "*.html"
    path "*.zip"
    script:
    """
    fastqc ${reads.join(' ')} --threads ${task.cpus} -o ./
    """
}

workflow {
    // CRITICAL FIX: flat:false
    read_pairs = channel.fromFilePairs(params.reads, flat:false)
    FASTQC(read_pairs)
}

