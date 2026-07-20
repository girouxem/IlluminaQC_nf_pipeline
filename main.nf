#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.reads          = params.reads          ?: "data/*_{r1,R1,r2,R2}.fastq.gz"
params.outdir         = params.outdir         ?: "${baseDir}/results"
params.busco_lineage  = params.busco_lineage  ?: 'bacteria_odb10'
params.min_contig_len = params.min_contig_len ?: 500
// Typing databases (set on CLI or leave null to skip)
params.kraken_db      = params.kraken_db      ?: null
params.virulence_db   = params.virulence_db   ?: null
// Resource defaults (override on CLI: --global_cpus 8 --global_memory '32 GB'
//   or per-stage: --cpus_spades 16 --mem_spades '64 GB')
params.global_cpus    = params.global_cpus    ?: 2
params.global_memory  = params.global_memory  ?: '4 GB'
params.cpus_spades    = params.cpus_spades    ?: null
params.mem_spades     = params.mem_spades     ?: null
params.cpus_busco     = params.cpus_busco     ?: null
params.mem_busco      = params.mem_busco      ?: null
params.cpus_kraken2   = params.cpus_kraken2   ?: null
params.mem_kraken2    = params.mem_kraken2    ?: null
params.cpus_quast     = params.cpus_quast     ?: null
params.mem_quast      = params.mem_quast      ?: null


// Helper: decide organism from Kraken2 report
def decideOrganism(String kraken_report_path) {
    def lines = new File(kraken_report_path).readLines()
    def species = lines.find { it.contains("\tS\t") } ?: ''
    if (species.contains("Mycobacterium")) return "mbovis"
    if (species.contains("Salmonella"))    return "salmonella"
    if (species.contains("Listeria"))      return "listeria"
    if (species.toLowerCase().contains("virus")) return "virus"
    return "bacteria"
}


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
        echo ${reads} | tr ' ' '\\n' >> paircheck_${sample_id}.txt
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

process RETENTION_TRACK {
    tag "$sample_id"
    publishDir "${params.outdir}/retention", mode: 'copy'
    cpus   { params.global_cpus }
    memory { params.global_memory }
    input:
    tuple val(sample_id), path(raw_pair), path(trim_R1), path(trim_R2)
    output:
    path "retention_${sample_id}.csv"
    script:
    """
    #!/bin/bash
    set -euo pipefail
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


// ---------- 7. SPAdes Assembly ----------
process SPADES {
    errorStrategy 'ignore'
    maxRetries 0
    publishDir "${params.outdir}/assembly", mode: 'copy'
    tag "$sample_id"
    cpus   { params.cpus_spades ?: params.global_cpus }
    memory { params.mem_spades  ?: params.global_memory }
    input:
    tuple val(sample_id), path(R1), path(R2)
    output:
    tuple val(sample_id), path("${sample_id}_contigs.fasta"), optional: true, emit: contigs
    path "${sample_id}_assembly_*txt", optional: true
    script:
    """
    #!/bin/bash
    set -euo pipefail
    if [[ ! -s "${R1}" || ! -s "${R2}" ]]; then
        echo "No trimmed data for ${sample_id}" > ${sample_id}_assembly_skipped.txt
        exit 0
    fi
    spades.py \\
        -1 "${R1}" \\
        -2 "${R2}" \\
        --isolate \\
        -o "${sample_id}_spades" \\
        --threads ${task.cpus} \\
        --memory ${(task.memory.toMega()/1024) as int} \\
        2>&1 | tee spades_run_${sample_id}.log
    if [[ -f "${sample_id}_spades/contigs.fasta" ]]; then
        cp "${sample_id}_spades/contigs.fasta" "${sample_id}_contigs.fasta"
    else
        echo "No contigs generated for ${sample_id}" > ${sample_id}_assembly_failed.txt
    fi
    """
}


// ---------- Pre-filter contigs before BUSCO ----------
process FILTER_CONTIGS {
    tag "$sample_id"
    publishDir "${params.outdir}/assembly_filtered", mode: 'copy'
    cpus   { params.global_cpus }
    memory { params.global_memory }

    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), path("${sample_id}_contigs_filtered.fasta"), emit: filtered_contigs

    script:
    """
    awk -v MIN=${params.min_contig_len} '
    BEGIN {RS=">"; ORS=""}
    NR>1 {
        idx = index(\$0, "\\n")
        if (idx > 0) {
            header = substr(\$0, 1, idx-1)
            seq = substr(\$0, idx+1)
            gsub(/\\n/, "", seq)
            if (length(seq) >= MIN) {
                print ">" header "\\n" seq "\\n"
            }
        }
    }
    ' ${contigs} > ${sample_id}_contigs_filtered.fasta
    """
}


// ---------- 8. QUAST Assembly QC ----------
process QUAST {
    tag "$sample_id"
    publishDir "${params.outdir}/quast", mode: 'copy'
    cpus   { params.cpus_quast ?: params.global_cpus }
    memory { params.mem_quast  ?: params.global_memory }
    input:
    tuple val(sample_id), path(contigs)
    output:
    path "${sample_id}_quast", emit: quast_reports
    script:
    """
    quast.py ${contigs} -o ${sample_id}_quast -t ${task.cpus}
    """
}


// ---------- 9. BUSCO Completeness ----------
process BUSCO {
    tag "$sample_id"
    publishDir "${params.outdir}/busco", mode: 'copy'
    errorStrategy 'ignore'
    maxRetries 0
    cpus   { params.cpus_busco ?: params.global_cpus }
    memory { params.mem_busco  ?: params.global_memory }
    input:
    tuple val(sample_id), path(contigs)
    output:
    path("${sample_id}_busco"), optional: true, emit: busco_reports
    path("${sample_id}_busco_failed.txt"), optional: true
    script:
    """
    #!/bin/bash
    set -euo pipefail
    if [[ ! -s "${contigs}" ]]; then
        echo "No contigs for ${sample_id}" > ${sample_id}_busco_failed.txt
        exit 0
    fi
    busco -i ${contigs} \\
          -o ${sample_id}_busco \\
          -l ${params.busco_lineage} \\
          -m genome \\
          --cpu ${task.cpus} \\
          2>&1 | tee busco_run_${sample_id}.log
    if [[ ! -d ${sample_id}_busco ]]; then
        echo "BUSCO failed for ${sample_id}" > ${sample_id}_busco_failed.txt
    fi
    """
}

// ---------- Collect BUSCO results into a combined table ----------
process BUSCO_SUMMARY {
    publishDir "${params.outdir}/busco", mode: 'copy'
    input:
    path busco_dirs
    output:
    path "busco_summary.tsv"
    script:
    """
    #!/bin/bash
    printf "sample\\tcomplete_pct\\tsingle_pct\\tduplicated_pct\\tfragmented_pct\\tmissing_pct\\ttotal_genes\\n" > busco_summary.tsv
    for d in ${busco_dirs}; do
        sample=\$(basename \$d | sed 's/_busco//')
        # Find the summary file - check common locations
        summary=""
        if [ -f "\$d/short_summary.txt" ]; then
            summary="\$d/short_summary.txt"
        fi
        if [ -z "\$summary" ]; then
            summary=\$(ls \$d/run_*/short_summary*.txt 2>/dev/null | head -1)
        fi
        if [ -z "\$summary" ]; then
            summary=\$(ls \$d/*/short_summary*.txt 2>/dev/null | head -1)
        fi
        if [ -z "\$summary" ]; then
            summary=\$(ls \$d/short_summary*.txt 2>/dev/null | head -1)
        fi
        if [ -z "\$summary" ]; then
            printf "%s\\tNO_SUMMARY_FILE\\tNA\\tNA\\tNA\\tNA\\tNA\\n" "\$sample" >> busco_summary.tsv
            continue
        fi
        # Match: C:19.3%[S:18.4%,D:0.9%],F:3.6%,M:77.1%,n:440
        line=\$(grep -E 'C:[0-9]' "\$summary" 2>/dev/null | head -1)
        if [ -z "\$line" ]; then
            printf "%s\\tNO_DATA_LINE\\tNA\\tNA\\tNA\\tNA\\tNA\\n" "\$sample" >> busco_summary.tsv
            continue
        fi
        complete=\$(echo "\$line" | sed -n 's/.*C:\\([0-9.]*\\)%.*/\\1/p')
        single=\$(echo "\$line"  | sed -n 's/.*S:\\([0-9.]*\\)%.*/\\1/p')
        dup=\$(echo "\$line"     | sed -n 's/.*D:\\([0-9.]*\\)%.*/\\1/p')
        frag=\$(echo "\$line"    | sed -n 's/.*F:\\([0-9.]*\\)%.*/\\1/p')
        miss=\$(echo "\$line"    | sed -n 's/.*M:\\([0-9.]*\\)%.*/\\1/p')
        total=\$(echo "\$line"   | sed -n 's/.*n:\\([0-9]*\\).*/\\1/p')
        printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" "\$sample" "\${complete:-PARSE_FAIL}" "\${single:-PARSE_FAIL}" "\${dup:-PARSE_FAIL}" "\${frag:-PARSE_FAIL}" "\${miss:-PARSE_FAIL}" "\${total:-PARSE_FAIL}" >> busco_summary.tsv
    done
    """
}

// ---------- Kraken2 classification ----------
process KRAKEN2 {
    tag "$sample_id"
    publishDir "${params.outdir}/kraken2", mode: 'copy'
    errorStrategy 'ignore'
    maxRetries 0
    cpus   { params.cpus_kraken2 ?: params.global_cpus }
    memory { params.mem_kraken2  ?: params.global_memory }
    input:
    tuple val(sample_id), path(contigs)
    path kraken_db                         // Captures the staged database directory
    output:
    tuple val(sample_id), path("${sample_id}_kraken2.report"), optional: true, emit: kraken_reports
    path "${sample_id}_kraken2.out", optional: true, emit: kraken_logs

    when:
    params.kraken_db != null
    script:
    """
    #!/bin/bash
    set -euo pipefail
    kraken2 --db ${params.kraken_db} \\
            --report ${sample_id}_kraken2.report \\
            ${contigs} > ${sample_id}_kraken2.out 2>&1
    """
}
// ---------- TB-Profiler (Mycobacterium bovis) ----------
process TB_PROFILER {
    tag "$sample_id"
    publishDir "${params.outdir}/tbprofiler", mode: 'copy'
    errorStrategy = 'ignore'
    maxRetries = 0

    input:
    tuple val(sample_id), path(R1), path(R2)

    output:
    tuple val(sample_id), path("results/${sample_id}.results.json"), optional: true, emit: tbprofiler_json
    tuple val(sample_id), path("results/${sample_id}.results.txt"), optional: true, emit: tbprofiler_txt
    path "${sample_id}_tbprofiler_failed.txt", optional: true

    script:
    """
    #!/bin/bash
    set -euo pipefail

    tb-profiler profile \
        --read1 ${R1} \
        --read2 ${R2} \
        --prefix ${sample_id} \
        --db ${HOME}/databases/tbprofiler/tbdb \
        --txt \
        --threads ${task.cpus} 2>&1 | tee tbprofiler_run_${sample_id}.log || true

    if [[ ! -f results/${sample_id}.results.json ]]; then
        echo "TB-Profiler failed for ${sample_id}" > ${sample_id}_tbprofiler_failed.txt
    fi
    """
}


// ---------- SISTR2 (Salmonella) ----------
process SISTR2 {
    tag "$sample_id"
    publishDir "${params.outdir}/sistr2", mode: 'copy'
    errorStrategy 'ignore'
    maxRetries 0
    cpus   { params.global_cpus }
    memory { params.global_memory }
    input:
    tuple val(sample_id), path(contigs)
    output:
    tuple val(sample_id), path("${sample_id}_sistr2.json"), optional: true, emit: sistr2_json
    script:
    """
    sistr --qc -i ${contigs} -f json -o ${sample_id}_sistr2.json
    """
}
// ---------- AMRFinderPlus (Salmonella AMR) ----------
process AMRFINDER {
    tag "$sample_id"
    publishDir "${params.outdir}/amrfinder", mode: 'copy'
    errorStrategy 'ignore'
    maxRetries 0
    cpus   { params.global_cpus }
    memory { params.global_memory }
    input:
    tuple val(sample_id), path(contigs)
    output:
    tuple val(sample_id), path("${sample_id}_amrfinder.txt"), optional: true, emit: amr_txt
    script:
    """
    amrfinder -n ${contigs} -o ${sample_id}_amrfinder.txt
    """
}
// ---------- MLST (Listeria) ----------
process MLST {
    tag "$sample_id"
    publishDir "${params.outdir}/mlst", mode: 'copy'
    errorStrategy 'ignore'
    maxRetries 0
    cpus   { params.global_cpus }
    memory { params.global_memory }
    input:
    tuple val(sample_id), path(contigs)
    output:
    tuple val(sample_id), path("${sample_id}_mlst.txt"), optional: true, emit: mlst_txt
    script:
    """
    mlst ${contigs} > ${sample_id}_mlst.txt
    """
}
// ---------- VirulenceFinder (Listeria) ----------
process VIRULENCEFINDER {
    tag "$sample_id"
    publishDir "${params.outdir}/virulencefinder", mode: 'copy'
    errorStrategy 'ignore'
    maxRetries 0
    cpus   { params.global_cpus }
    memory { params.global_memory }
    input:
    tuple val(sample_id), path(contigs)
    path virulence_db                     // Receives the database via an auto-mounting Nextflow path channel
    output:
    tuple val(sample_id), path("${sample_id}_virulence.json"), optional: true, emit: vir_json
    when:
    params.virulence_db != null
    script:
    """
    #!/bin/bash
    set -euo pipefail
    virulencefinder.py \\
       -i ${contigs} \\
       -o ${sample_id}_vf \\
       -p ${params.virulence_db} \\
       -x \\
       --json
    cp ${sample_id}_vf/*.json ${sample_id}_virulence.json 2>/dev/null || touch ${sample_id}_virulence.json
    """
}

// ---------- MULTIQC (now includes QUAST + BUSCO) ----------
process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    cpus   { params.global_cpus }
    memory { params.global_memory }
    input:
    path all_reports
    output:
    path "multiqc_report.html"
    script:
    """
    multiqc . -o ./
    """
}

/*
 * NORMAL FULL WORKFLOW
 * Run with:
 *   nextflow run .
 */

// Full pipeline
workflow {
    log.info "Output dir       : ${params.outdir}"
    log.info "Reads glob       : ${params.reads}"
    log.info "BUSCO lineage    : ${params.busco_lineage}"
    log.info "Min contig len   : ${params.min_contig_len}"
    log.info "Global CPUs      : ${params.global_cpus}"
    log.info "Global memory    : ${params.global_memory}"

    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    SANITY_CHECK(raw_pairs)
    raw_fastqc = FASTQC_RAW(raw_pairs)
    trimmed    = FASTP(raw_pairs)
    // trimmed.trimmed_reads is already (sample_id, R1, R2) — feed directly
    trimmed_fastqc = FASTQC_TRIMMED(trimmed.trimmed_reads)
    raw_reports     = raw_fastqc.fastqc_zip.map { id, f -> f }
    trimmed_reports = trimmed.fastp_json.map { id, f -> f }
    post_reports    = trimmed_fastqc.fastqc_zip_trim.map { id, f -> f }
    retention_input = raw_pairs.join(trimmed.trimmed_reads)
    RETENTION_TRACK(retention_input)
    // Assembly -> filter -> QUAST + BUSCO
    assemblies = SPADES(trimmed.trimmed_reads)
    filtered   = FILTER_CONTIGS(assemblies.contigs)
    quast_out  = QUAST(assemblies.contigs)
    busco_out  = BUSCO(filtered.filtered_contigs)
    // Collect all BUSCO dirs into one channel, then summarize
    busco_collected = busco_out.busco_reports.collect()
    BUSCO_SUMMARY(busco_collected)

    // Kraken2 classification for routing
    ch_kraken_db = params.kraken_db ? file(params.kraken_db, checkIfExists: true) : []
    kraken_out = KRAKEN2(assemblies.contigs, ch_kraken_db)
    // Route samples to typing tools based on Kraken2 result
    if (params.kraken_db) {
        // TB-Profiler: needs reads, run on M. bovis samples
        tb_inputs = trimmed.trimmed_reads.join(kraken_out.kraken_reports)
                     .filter { id, r1, r2, rep -> decideOrganism(rep.toString()) == 'mbovis' }
                     .map    { id, r1, r2, rep -> tuple(id, r1, r2) }
        TB_PROFILER(tb_inputs)
        // SISTR2 + AMRFinderPlus: run on Salmonella contigs
        sal_inputs = assemblies.contigs.join(kraken_out.kraken_reports)
                      .filter { id, c, rep -> decideOrganism(rep.toString()) == 'salmonella' }
                      .map    { id, c, rep -> tuple(id, c) }
        SISTR2(sal_inputs)
        AMRFINDER(sal_inputs)
        // MLST + VirulenceFinder: run on Listeria contigs
        lis_inputs = assemblies.contigs.join(kraken_out.kraken_reports)
                      .filter { id, c, rep -> decideOrganism(rep.toString()) == 'listeria' }
                      .map    { id, c, rep -> tuple(id, c) }
        MLST(lis_inputs)
        ch_virulence_db = params.virulence_db ? file(params.virulence_db, checkIfExists: true) : []
        VIRULENCEFINDER(lis_inputs, ch_virulence_db)

    }
    // MultiQC: now includes QUAST + BUSCO outputs
    quast_reports = quast_out.quast_reports.map { p -> p }
    busco_reports = busco_out.busco_reports.map { p -> p }
    all_reports = raw_reports
        .mix(trimmed_reports, post_reports, quast_reports, busco_reports)
        .collect()
    MULTIQC(all_reports)
}

/*
 * TEST 1: sanity check only
 * Run with:
 *   nextflow run . -entry test_SANITY_CHECK
 */
workflow test_SANITY_CHECK {
    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    raw_pairs
        .map { it.toString() }
        .view()
    SANITY_CHECK(raw_pairs)
}

/*
 * TEST 2: raw FastQC only
 * Run with:
 *   nextflow run . -entry test_FASTQC_RAW
 */
workflow test_FASTQC_RAW {
    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    SANITY_CHECK(raw_pairs)
    raw_pairs
        .map { it.toString() }
        .view()
    FASTQC_RAW(raw_pairs)
}

/*
 * TEST 3: fastp only
 * Run with:
 *   nextflow run . -entry test_FASTP
 */
workflow test_FASTP {
    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    SANITY_CHECK(raw_pairs)
    raw_pairs
        .map { it.toString() }
        .view()
    trimmed = FASTP(raw_pairs)
    trimmed.trimmed_reads
        .map { it.toString() }
        .view()
}

/*
 * TEST 4: FastQC on trimmed reads only
 * Run with:
 *   nextflow run . -entry test_FASTQC_TRIMMED
 */
workflow test_FASTQC_TRIMMED {
    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    trimmed = FASTP(raw_pairs)
    // trimmed.trimmed_reads is already (sample_id, R1, R2) — feed directly
    trimmed.trimmed_reads
        .map { it.toString() }
        .view()
    FASTQC_TRIMMED(trimmed.trimmed_reads)
}

/*
 * TEST 5: retention only
 * Run with:
 *   nextflow run . -entry test_RETENTION_TRACK
 */
workflow test_RETENTION_TRACK {
    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    trimmed = FASTP(raw_pairs)
    retention_input = raw_pairs.join(trimmed.trimmed_reads)
    retention_input
        .map { it.toString() }
        .view()
    RETENTION_TRACK(retention_input)
}

/*
 * TEST 6: SPAdes only
 * Run with:
 *   nextflow run . -entry test_SPADES
 */
workflow test_SPADES {
    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    trimmed = FASTP(raw_pairs)
    trimmed.trimmed_reads
        .map { it.toString() }
        .view()
    SPADES(trimmed.trimmed_reads)
}

/*
 * TEST 7: QUAST only
 * Run with:
 *   nextflow run . -entry test_QUAST
 */
workflow test_QUAST {
    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    trimmed = FASTP(raw_pairs)
    assemblies = SPADES(trimmed.trimmed_reads)
    assemblies.contigs
        .map { it.toString() }
        .view()
    QUAST(assemblies.contigs)
}

/*
 * TEST 8: BUSCO only
 * Run with:
 *   nextflow run . -entry test_BUSCO --busco_lineage enterobacterales_odb10
 */
workflow test_BUSCO {
    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    trimmed = FASTP(raw_pairs)
    assemblies = SPADES(trimmed.trimmed_reads)
    filtered = FILTER_CONTIGS(assemblies.contigs)
    busco_out = BUSCO(filtered.filtered_contigs)
    busco_collected = busco_out.busco_reports.collect()
    BUSCO_SUMMARY(busco_collected)
}

workflow test_BUSCO_SUMMARY{
    def real_busco_ch = channel.fromPath("${params.outdir}/busco/*_busco", type: 'dir', checkIfExists: true)
        .collect()
    // Crucial: gathers all folders into a single list
    BUSCO_SUMMARY(real_busco_ch)
}
/*
 * TEST 9: KRAKEN2 only
 * Run with:
 *   nextflow run . -entry test_KRAKEN2 --kraken_db ~/databases/kraken2
 */
workflow test_KRAKEN2 {
    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    trimmed = FASTP(raw_pairs)
    assemblies = SPADES(trimmed.trimmed_reads)
    ch_kraken_db = file(params.kraken_db, checkIfExists: true)
    assemblies.contigs
        .map { it.toString() }
        .view()
    KRAKEN2(assemblies.contigs, ch_kraken_db)
}

/*
 * TEST 10: TB-Profiler only (reuses cached FASTP results)
 * Run with:
 *   nextflow run . -entry test_TB_PROFILER -resume
 */
workflow test_TB_PROFILER {

    raw_pairs = channel.fromFilePairs(params.reads, flat: false)

    trimmed = FASTP(raw_pairs)

    // Show what TB-Profiler will receive
    trimmed.trimmed_reads
        .map { it.toString() }
        .view()

    TB_PROFILER(trimmed.trimmed_reads)
}
