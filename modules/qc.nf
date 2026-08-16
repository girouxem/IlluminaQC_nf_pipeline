// modules/qc.nf
// Contains: SANITY_CHECK, FASTQC_RAW, FASTP, FASTQC_TRIMMED

process SANITY_CHECK {
    publishDir "${params.outdir}/logs", mode: 'copy'
    cpus   { params.global_cpus }
    memory { params.global_memory }
    input:
    tuple val(sample_id), path(reads)
    output:
    path "paircheck_${sample_id}.txt"
    script:
    """
    echo "Sample ID: ${sample_id}" > paircheck_${sample_id}.txt
    read_count=\$(echo ${reads} | wc -w)
    if [ "\$read_count" -eq 2 ]; then
        echo "OK ${sample_id} has two files:" >> paircheck_${sample_id}.txt
        echo ${reads} | tr ' ' '\\n' >> paircheck_${sample_id}.txt
    else
        echo "ERROR ${sample_id} has \$read_count files (expected 2)" >> paircheck_${sample_id}.txt
        exit 1
    fi
    """
}

process FASTQC_RAW {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc_raw", mode: 'copy'
    cpus   { params.global_cpus }
    memory { params.global_memory }
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
    cpus   { params.global_cpus }
    memory { params.global_memory }
    input:
    tuple val(sample_id), path(reads)
    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R1.fastq.gz"),
                          path("${sample_id}_trimmed_R2.fastq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("${sample_id}_fastp_report.html"), emit: fastp_html
    tuple val(sample_id), path("${sample_id}_fastp.json"), emit: fastp_json
    script:
    """
    fastp \\
      -i ${reads[0]} \\
      -I ${reads[1]} \\
      -o ${sample_id}_trimmed_R1.fastq.gz \\
      -O ${sample_id}_trimmed_R2.fastq.gz \\
      --json ${sample_id}_fastp.json \\
      --html ${sample_id}_fastp_report.html \\
      --thread ${task.cpus}
    """
}

process FASTQC_TRIMMED {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc_trimmed", mode: 'copy'
    cpus   { params.global_cpus }
    memory { params.global_memory }
    input:
    tuple val(sample_id), path(read1), path(read2)
    output:
    tuple val(sample_id), path("*.zip"), emit: fastqc_zip_trim
    tuple val(sample_id), path("*.html"), emit: fastqc_html_trim
    script:
    """
    fastqc ${read1} ${read2} --threads ${task.cpus} -o ./
    """
}
