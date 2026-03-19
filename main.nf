#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
* Illumina QC pipeline:
* 1. FASTQC (raw)
* 2. fastq trimming
* 3. FASTQC (trimmed)
* 4. MultiQC aggregation
*/
params.reads  = "data/*_{r1,R1,r2,R2}.fastq.gz"
params.outdir = "${baseDir}/results"

process SANITY_CHECK {

    publishDir "${params.outdir}/logs", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "paircheck_${sample_id}.txt"

    script:
    """
    echo "Sample ID: ${sample_id}" > paircheck_${sample_id}.txt
    read_count=\$(echo ${reads} | wc -w)
    if [ "\$read_count" -eq 2 ]; then
        echo "✅  ${sample_id} has two files:" >> paircheck_${sample_id}.txt
        echo ${reads} | tr ' ' '\\n' >> paircheck_${sample_id}.txt
    else
        echo "❌  ${sample_id} has \$read_count files (expected 2)" >> paircheck_${sample_id}.txt
        echo ${reads} | tr ' ' '\\n' >> paircheck_${sample_id}.txt
        exit 1
    fi
    """
}



process FASTQC_RAW {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc_raw", mode: 'copy'
    input:
    tuple val(sample_id), path(reads)
    output:
    tuple val(sample_id), path("*.zip"), emit: fastqc_zip
    tuple val(sample_id), path("*.html"), emit: fastqc_html
    script:
    """
    fastqc ${reads.join(' ')} --threads ${task.cpus} -o ./
    """
}
process FASTP {
    tag "$sample_id"
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    input:
    tuple val(sample_id), path(reads)
    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R1.fastq.gz"), 
                          path("${sample_id}_trimmed_R2.fastq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("${sample_id}_fastp_report.html"), emit: fastp_html
    tuple val(sample_id), path("${sample_id}_fastp.json"), emit: fastp_json
    script:
    """
    fastp \
      -i ${reads[0]} \
      -I ${reads[1]} \
      -o ${sample_id}_trimmed_R1.fastq.gz \
      -O ${sample_id}_trimmed_R2.fastq.gz \
      --json ${sample_id}_fastp.json \
      --html ${sample_id}_fastp_report.html \
      --thread ${task.cpus}
    """
}
process FASTQC_TRIMMED {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc_trimmed", mode: 'copy'
    input:
    tuple val(sample_id), path(read1), path(read2)
    output:
    tuple val(sample_id), path("*.zip"), emit: fastqc_zip_trim
    tuple val(sample_id), path("*.html"), emit: fastqc_html_trim
    script:
    """
    fastqc ${read1} ==threads ${task.cpus} -o ./
    """
}

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    input:
    path all_reports
    output:
    path "multiqc_report.html"
    script:
    """
    multiqc . -o ./
    """
}


workflow {
    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    SANITY_CHECK(raw_pairs)
    raw_fastqc     = FASTQC_RAW(raw_pairs)
    trimmed        = FASTP(raw_pairs)
    trimmed_fastqc = FASTQC_TRIMMED(trimmed.trimmed_reads)
    // Take only the file paths (ignore sample IDs)
    raw_reports     = raw_fastqc.fastqc_zip.map { id, f -> f }
    trimmed_reports = trimmed.fastp_json.map   { id, f -> f }
    post_reports    = trimmed_fastqc.fastqc_zip_trim.map { id, f -> f }
    // Merge all report files, then collect them
    all_reports = raw_reports
        .mix(trimmed_reports, post_reports)
        .collect()
    MULTIQC(all_reports)
}