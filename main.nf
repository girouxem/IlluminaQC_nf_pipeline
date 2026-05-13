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

process RETENTION_TRACK {

    tag "$sample_id"
    publishDir "${params.outdir}/retention", mode: 'copy'

    /*
     * Input structure from the join():
     * [ sample_id , [raw_R1, raw_R2] , trim_R1 , trim_R2 ]
     */
    input:
    tuple val(sample_id), path(raw_pair), path(trim_R1), path(trim_R2)

    output:
    path "retention_${sample_id}.csv"

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # choose the first mate from the raw pair
    raw_R1=\$(echo ${raw_pair} | cut -d' ' -f1)

    raw_count=\$(zgrep -c '^@' "\${raw_R1}")
    trim_count=\$(zgrep -c '^@' "${trim_R1}")

    raw_bases=\$(zcat "\${raw_R1}" | paste - - - - | cut -f2 | tr -d '\\n' | wc -c)
    trim_bases=\$(zcat "${trim_R1}" | paste - - - - | cut -f2 | tr -d '\\n' | wc -c)

    {
      echo "sample,reads_before,reads_after,bases_before,bases_after"
      echo "${sample_id},\${raw_count},\${trim_count},\${raw_bases},\${trim_bases}"
    } > "retention_${sample_id}.csv"
    """
}


// ---------- 7.  SPAdes Assembly ----------
process SPADES {

    // directives (in any order among themselves)
    errorStrategy = 'ignore'
    maxRetries    = 0
    cpus          = 8
    memory        = '16 GB'
    publishDir "${params.outdir}/assembly", mode: 'copy'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(R1), path(R2)

    output:
    tuple val(sample_id), path("contigs.fasta"), emit: contigs

    script:
    """
    #!/bin/bash
    set -euo pipefail

    if [[ ! -s "${R1}" || ! -s "${R2}" ]]; then
        echo "No trimmed data for ${sample_id}" > ${sample_id}_assembly_skipped.txt
        exit 0
    fi

    spades.py \
        -1 "${R1}" \
        -2 "${R2}" \
        --isolate \
        -o "${sample_id}_spades" \
        --threads ${task.cpus} \
        --memory ${task.memory.toMega()/1024} \
        2>&1 | tee spades_run_${sample_id}.log || true

    if [[ -f "${sample_id}_spades/contigs.fasta" ]]; then
        cp "${sample_id}_spades/contigs.fasta" .
    else
        echo "No contigs generated for ${sample_id}" > ${sample_id}_assembly_failed.txt
    fi
    """
}



// ---------- 8.  QUAST Assembly QC ----------
process QUAST {

    tag "$sample_id"
    publishDir "${params.outdir}/quast", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs)

    output:
    path "${sample_id}_quast", emit: quast_reports

    script:
    """
    quast.py ${contigs} -o ${sample_id}_quast -t ${task.cpus}
    """
}


// ---------- 9.  BUSCO Completeness ----------
// ---------- 9. BUSCO Completeness ----------
process BUSCO {

    tag "$sample_id"
    publishDir "${params.outdir}/busco", mode: 'copy'
    errorStrategy = 'ignore'
    maxRetries = 0

    input:
    tuple val(sample_id), path(contigs)

    output:
    path("${sample_id}_busco") , optional: true , emit: busco_reports
    path("${sample_id}_busco_failed.txt") , optional: true

    script:
    """
    #!/bin/bash
    set -euo pipefail

    if [[ ! -s "${contigs}" ]]; then
        echo "No contigs for ${sample_id}" > ${sample_id}_busco_failed.txt
        exit 0
    fi

    busco -i ${contigs} \
          -o ${sample_id}_busco \
          -l ${params.busco_lineage} \
          -m genome \
          --cpu ${task.cpus} \
          2>&1 | tee busco_run_${sample_id}.log || true

    if [[ ! -d ${sample_id}_busco ]]; then
        echo "BUSCO failed for ${sample_id}" > ${sample_id}_busco_failed.txt
    fi
    """
}





workflow {
    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    SANITY_CHECK(raw_pairs)
    raw_fastqc     = FASTQC_RAW(raw_pairs)
    trimmed        = FASTP(raw_pairs)
    trimmed_fastqc = FASTQC_TRIMMED(trimmed.trimmed_reads)
    // Prepare channel: link each sample’s raw+trimmed reads
    // Take only the file paths (ignore sample IDs)
    raw_reports     = raw_fastqc.fastqc_zip.map { id, f -> f }
    trimmed_reports = trimmed.fastp_json.map   { id, f -> f }
    post_reports    = trimmed_fastqc.fastqc_zip_trim.map { id, f -> f }
    // Merge all report files, then collect them
    all_reports = raw_reports
        .mix(trimmed_reports, post_reports)
        .collect()
    retention_input = raw_pairs.join(trimmed.trimmed_reads)
    RETENTION_TRACK(retention_input)

    // ---------- Assembly workflow ----------
    trimmed_reads_flat = trimmed.trimmed_reads  // emits [id, R1, R2]
    assemblies = SPADES(trimmed_reads_flat)
    quast_out  = QUAST(assemblies.contigs)
    busco_out  = BUSCO(assemblies.contigs)

    MULTIQC(all_reports)
}