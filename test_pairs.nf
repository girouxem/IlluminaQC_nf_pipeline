#!/usr/bin/env nextflow
nextflow.enable.dsl=2
workflow {
    channel
        .fromFilePairs("data/*_{r1,r2}.fastq.gz", flat: true)
        .view()
}

