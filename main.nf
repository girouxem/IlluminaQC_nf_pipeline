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
// ---------- Kraken2 classification (routes samples to typing tools) ----------
process KRAKEN2 {
    tag "$sample_id"
    publishDir "${params.outdir}/kraken2", mode: 'copy'
    errorStrategy 'ignore'
    maxRetries 0

    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), path("${sample_id}_kraken2.report"), optional: true, emit: kraken_reports

    when:
    params.kraken_db != null

    script:
    """
    #!/bin/bash
    set -euo pipefail

    kraken2 --db ${params.kraken_db} \
            --memory-mapping \
            --threads ${task.cpus} \
            --report ${sample_id}_kraken2.report \
            ${contigs} > ${sample_id}_kraken2.out 2>&1

    if [[ ! -f ${sample_id}_kraken2.report ]]; then
        echo "Kraken2 failed for ${sample_id}" > ${sample_id}_kraken2.report
    fi
    """

}
// ---------- TB-Profiler (Mycobacterium bovis) ----------
process TB_PROFILER {
    tag "$sample_id"
    publishDir "${params.outdir}/tbprofiler", mode: 'copy'
    errorStrategy 'ignore'
    maxRetries 0

    input:
    tuple val(sample_id), path(R1), path(R2)

    output:
    tuple val(sample_id), path("${sample_id}.results.json"), optional: true, emit: tbprofiler_json
    tuple val(sample_id), path("${sample_id}.results.txt"), optional: true, emit: tbprofiler_txt
    tuple val(sample_id), path("${sample_id}_tbprofiler_mqc.tsv"), optional: true, emit: tbprofiler_mqc
    path "${sample_id}_tbprofiler_failed.txt", optional: true

        script:
    """
    #!/bin/bash

    tb-profiler profile \
        --read1 ${R1} \
        --read2 ${R2} \
        --prefix ${sample_id} \
        --db ${HOME}/databases/tbprofiler/tbdb \
        --txt \
        --threads ${task.cpus} 2>&1 | tee tbprofiler_run_${sample_id}.log || true

    LINEAGE="NA"
    SUBLINEAGE="NA"
    DRUGS="None"

    if [[ -f results/${sample_id}.results.json ]]; then
        cp results/${sample_id}.results.json . 2>/dev/null || true
        cp results/${sample_id}.results.txt . 2>/dev/null || true

                # Parse v6.x JSON with python3 -c (no heredoc, no file)
        python3 -c "
import json
d=json.load(open('${sample_id}.results.json'))
lin=d.get('lineage',[])
lineage=lin[0].get('lineage','NA') if lin else 'NA'
sublineage=d.get('sub_lineage','NA')
if not sublineage or sublineage=='NA':
    sublineage=lin[0].get('family','NA') if lin else 'NA'
dr=d.get('dr_variants',[])
drugs=[]
for v in dr:
    for ann in v.get('annotation',[]):
        if ann.get('drug'): drugs.append(ann['drug'])
drugs=sorted(set(drugs))
drugs_str=';'.join(drugs) if drugs else 'None'
with open('${sample_id}_tbprofiler_mqc.tsv','w') as out:
    out.write('Sample\\tLineage\\tSublineage\\tDrug_Resistance\\n')
    out.write('${sample_id}\\t'+lineage+'\\t'+sublineage+'\\t'+drugs_str+'\\n')
" 2>/dev/null || true


        if [[ ! -f ${sample_id}_tbprofiler_mqc.tsv ]]; then
            printf "Sample\\tLineage\\tSublineage\\tDrug_Resistance\\n" > ${sample_id}_tbprofiler_mqc.tsv
            printf "${sample_id}\\tPARSE_FAIL\\tPARSE_FAIL\\tPARSE_FAIL\\n" >> ${sample_id}_tbprofiler_mqc.tsv
        fi
    else
        echo "TB-Profiler failed for ${sample_id}" > ${sample_id}_tbprofiler_failed.txt
        printf "Sample\\tLineage\\tSublineage\\tDrug_Resistance\\n" > ${sample_id}_tbprofiler_mqc.tsv
        printf "${sample_id}\\tNA\\tNA\\tNone\\n" >> ${sample_id}_tbprofiler_mqc.tsv
    fi
    """
}
// ---------- Combine TB-Profiler results into one table ----------
process TB_SUMMARY {
    publishDir "${params.outdir}/tbprofiler", mode: 'copy'
    errorStrategy 'ignore'

    input:
    path tb_tsvs

    output:
    path "tbprofiler_summary.tsv", optional: true, emit: tb_summary

    script:
    """
    #!/bin/bash
    # Always create the file, even if empty
    printf "Sample\\tLineage\\tSublineage\\tDrug_Resistance\\n" > tbprofiler_summary.tsv

    # tb_tsvs may be a single file or a list  handle both
    for f in ${tb_tsvs}; do
        if [ -f "\$f" ]; then
            # Skip header row, append data
            tail -n +2 "\$f" >> tbprofiler_summary.tsv
        fi
    done
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
    tuple val(sample_id), path("${sample_id}_sistr2.tsv"), optional: true, emit: sistr2_tsv

    script:
    """
    #!/bin/bash

    # -i needs BOTH fasta path AND genome name
    sistr -i ${contigs} ${sample_id} -o ${sample_id}_sistr2 -f json 2>&1 | tee sistr_run_${sample_id}.log || true

    # sistr writes to <output>_predictions.json  rename it
    if [ -f ${sample_id}_sistr2_predictions.json ]; then
        cp ${sample_id}_sistr2_predictions.json ${sample_id}_sistr2.json
    fi

        # Generate TSV from JSON
    if [ -f ${sample_id}_sistr2.json ]; then
        python3 -c "
import json
data=json.load(open('${sample_id}_sistr2.json'))
# SISTR returns a LIST  take first element
d = data[0] if isinstance(data, list) and data else {}
# Try common serovar key names
serovar = d.get('serovar', d.get('serovar_predicted', d.get('serovar_antigenic', 'NA')))
# If still NA, build from antigenic_formula
if serovar == 'NA' and d.get('antigenic_formula'):
    serovar = 'antigenic: ' + str(d['antigenic_formula'])
serogroup = d.get('serogroup', 'NA')
h1 = d.get('h1', 'NA')
h2 = d.get('h2', 'NA')
cgmlst = d.get('cgmlst_ST', 'NA')
subspecies = d.get('cgmlst_subspecies', 'NA')
with open('${sample_id}_sistr2.tsv','w') as out:
    out.write('Sample\\tSerovar\\tSerogroup\\tH1\\tH2\\tCGMLST_ST\\tSubspecies\\n')
    out.write('${sample_id}\\t'+str(serovar)+'\\t'+str(serogroup)+'\\t'+str(h1)+'\\t'+str(h2)+'\\t'+str(cgmlst)+'\\t'+str(subspecies)+'\\n')
" 2>/dev/null || true
    else
        echo "SISTR2 failed for ${sample_id}" > ${sample_id}_sistr2_failed.txt
    fi

    """
}

// ---------- Combine SISTR2 results into one table ----------
process SISTR_SUMMARY {
    publishDir "${params.outdir}/sistr2", mode: 'copy'
    errorStrategy = 'ignore'

    input:
    path sistr_tsvs

    output:
    path "sistr2_summary.tsv", optional: true, emit: sistr_summary

    script:
    """
    #!/bin/bash
    printf "Sample\\tSerovar\\tSerogroup\\tH1\\tH2\\tCGMLST_ST\\tSubspecies\\n" > sistr2_summary.tsv
    for f in ${sistr_tsvs}; do
        if [ -f "\$f" ]; then
            tail -n +2 "\$f" >> sistr2_summary.tsv
        fi
    done
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
    tuple val(sample_id), path("${sample_id}_amrfinder.tsv"), optional: true, emit: amr_tsv
    script:
    """
    amrfinder -n ${contigs} -o ${sample_id}_amrfinder.tsv || true
    if [[ ! -s ${sample_id}_amrfinder.tsv ]]; then
        echo "Protein identifier,Contig id,Start,End,Strand,Gene symbol,Sequence name,Method,Scope,Element type,Element subtype,Class,Subclass,Method,Target length,Reference length,Coverage,Identity,Match type,PMID,Alignment length" > ${sample_id}_amrfinder.tsv
    fi
    """
}
// ---------- Combine AMRFinderPlus results into one table ----------
process AMR_SUMMARY {
    publishDir "${params.outdir}/amrfinder", mode: 'copy'
    errorStrategy 'ignore'

    input:
    path amr_tsvs

    output:
    path "amrfinder_summary.tsv", optional: true, emit: amr_summary

        script:
    """
    #!/bin/bash
    printf "Sample\\tAMR_Gene_Count\\tDrug_Classes\\tGene_Symbols\\n" > amrfinder_summary.tsv

    for f in ${amr_tsvs}; do
        if [ -f "\$f" ]; then
            # Get sample name using bash string replacement (avoids sed)
            f_name=\$(basename "\$f")
            sample="\${f_name%_amrfinder.tsv}"
            
            # Count data rows (skip header)
            count=\$(tail -n +2 "\$f" | wc -l)
            
            # Extract drug classes and gene symbols, join with ';' (no trailing semicolon)
            drugs=\$(tail -n +2 "\$f" | cut -f12 | sort -u | paste -sd ';' -)
            genes=\$(tail -n +2 "\$f" | cut -f6 | sort -u | paste -sd ';' -)
            
            printf "%s\\t%s\\t%s\\t%s\\n" "\$sample" "\$count" "\${drugs:-None}" "\${genes:-None}" >> amrfinder_summary.tsv
        fi
    done
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
    tuple val(sample_id), path("${sample_id}_mlst.tsv"), optional: true, emit: mlst_tsv
    script:
    """
    mlst ${contigs} > ${sample_id}_mlst.tsv || true
    if [[ ! -s ${sample_id}_mlst.tsv ]]; then
        echo "No MLST scheme found" > ${sample_id}_mlst.tsv
    fi
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

    output:
    tuple val(sample_id), path("${sample_id}_virulence.json"), optional: true, emit: vir_json
    tuple val(sample_id), path("${sample_id}_virulence.txt"), optional: true, emit: vir_txt

    when:
    params.virulence_db != null

    script:
    """
    #!/bin/bash
    set -euo pipefail

    if [[ ! -s "${contigs}" ]]; then
        echo "No contigs for ${sample_id}" > ${sample_id}_virulence.txt
        echo "{}" > ${sample_id}_virulence.json
        exit 0
    fi

    virulencefinder.py \
       -i ${contigs} \
       -o ${sample_id}_vf \
       -p ${params.virulence_db} \
       -x \
       --json 2>&1 | tee virulencefinder_run_${sample_id}.log || true

    cp ${sample_id}_vf/*.json ${sample_id}_virulence.json 2>/dev/null || echo "{}" > ${sample_id}_virulence.json
    cp ${sample_id}_vf/*.txt ${sample_id}_virulence.txt 2>/dev/null || echo "No virulence genes found" > ${sample_id}_virulence.txt
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
    multiqc . -o ./ -c ${baseDir}/multiqc_config.yaml
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
    kraken_out = KRAKEN2(assemblies.contigs)
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
        // Collect SISTR2 per-sample TSVs and merge
        sistr_tsv       = SISTR2.out.sistr2_tsv.map { id, f -> f }.collect().ifEmpty([])
        sistr_summary_out = SISTR_SUMMARY(sistr_tsv)

        AMRFINDER(sal_inputs)
        // Collect AMRFinderPlus per-sample TSVs and merge
        amr_reports     = AMRFINDER.out.amr_tsv.map { id, f -> f }.collect().ifEmpty([])
        amr_summary_out = AMR_SUMMARY(amr_reports)

        // MLST + VirulenceFinder: run on Listeria contigs
        lis_inputs = assemblies.contigs.join(kraken_out.kraken_reports)
                      .filter { id, c, rep -> decideOrganism(rep.toString()) == 'listeria' }
                      .map    { id, c, rep -> tuple(id, c) }
        MLST(lis_inputs)
        VIRULENCEFINDER(lis_inputs)


    }
    // MultiQC: now includes QUAST + BUSCO outputs
    quast_reports = quast_out.quast_reports.map { p -> p }
    busco_reports = busco_out.busco_reports.map { p -> p }
    
        // Collect QC reports
    raw_reports     = raw_fastqc.fastqc_zip.map { id, f -> f }
    trimmed_reports = trimmed.fastp_json.map { id, f -> f }
    post_reports    = trimmed_fastqc.fastqc_zip_trim.map { id, f -> f }
    quast_reports   = quast_out.quast_reports.map { p -> p }
    busco_reports   = busco_out.busco_reports.map { p -> p }

    // Collect typing reports (empty if tools didn't run)
    // Collect TB-Profiler per-sample TSVs and merge into one table
    tb_mqc_reports  = TB_PROFILER.out.tbprofiler_mqc.map { id, f -> f }.collect().ifEmpty([])
    tb_summary_out  = TB_SUMMARY(tb_mqc_reports)

    // Other typing reports (unchanged)
    tb_reports      = TB_PROFILER.out.tbprofiler_json.map { id, f -> f }.ifEmpty([])
    tb_txt_reports  = TB_PROFILER.out.tbprofiler_txt.map { id, f -> f }.ifEmpty([])
    sistr_reports   = SISTR2.out.sistr2_json.map { id, f -> f }.ifEmpty([])
    sistr_tsv       = SISTR2.out.sistr2_tsv.map { id, f -> f }.ifEmpty([])
    amr_reports     = AMRFINDER.out.amr_tsv.map { id, f -> f }.ifEmpty([])
    mlst_reports    = MLST.out.mlst_tsv.map { id, f -> f }.ifEmpty([])
    vir_reports     = VIRULENCEFINDER.out.vir_json.map { id, f -> f }.ifEmpty([])

    // Feed the combined summary to MultiQC (not the individual per-sample files)
    all_reports = raw_reports
        .mix(trimmed_reports, post_reports, quast_reports, busco_reports)
        .mix(tb_reports, tb_txt_reports, tb_summary_out.tb_summary,
             sistr_reports, sistr_summary_out.sistr_summary,
             amr_summary_out.amr_summary,
             mlst_reports, vir_reports)
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
    assemblies.contigs
        .map { it.toString() }
        .view()
    KRAKEN2(assemblies.contigs)

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

/*
 * TEST 11: Kraken2 routing (reuses cached SPADES assemblies)
 * Run with:
 *   nextflow run . -entry test_KRAKEN2_ROUTING --kraken_db ~/databases/kraken2 -resume
 *
 * Verifies:
 *   1. Kraken2 runs on each assembly
 *   2. decideOrganism() reads each report correctly
 *   3. Each sample is routed to the correct typing tool
 */
workflow test_KRAKEN2_ROUTING {

    raw_pairs = channel.fromFilePairs(params.reads, flat: false)

    trimmed = FASTP(raw_pairs)

    assemblies = SPADES(trimmed.trimmed_reads)

    // Run Kraken2 on each assembly
    kraken_out = KRAKEN2(assemblies.contigs)

    // Show what Kraken2 produced
    kraken_out.kraken_reports
        .map { id, rep -> "KRAKEN2: ${id} -> ${rep}" }
        .view()

    // Join assemblies with kraken reports so we can route
    // Shape after join: [sample_id, contigs, kraken_report]
    routed = assemblies.contigs.join(kraken_out.kraken_reports)

    // Apply the routing function and show the decision
    routed
        .map { id, contigs, rep ->
            def org = decideOrganism(rep.toString())
            "ROUTING: ${id} -> ${org}"
        }
        .view()

    // --- TB-Profiler: needs reads, route M. bovis samples ---
    tb_inputs = trimmed.trimmed_reads.join(kraken_out.kraken_reports)
                 .filter { id, r1, r2, rep -> decideOrganism(rep.toString()) == 'mbovis' }
                 .map    { id, r1, r2, rep -> tuple(id, r1, r2) }

    tb_inputs
        .map { id, r1, r2 -> "TB_PROFILER will run on: ${id}" }
        .view()

    TB_PROFILER(tb_inputs)

    // --- SISTR2 + AMRFinderPlus: route Salmonella samples ---
    sal_inputs = assemblies.contigs.join(kraken_out.kraken_reports)
                  .filter { id, c, rep -> decideOrganism(rep.toString()) == 'salmonella' }
                  .map    { id, c, rep -> tuple(id, c) }

    sal_inputs
        .map { id, c -> "SISTR2 + AMRFINDER will run on: ${id}" }
        .view()

    SISTR2(sal_inputs)
    AMRFINDER(sal_inputs)

    // --- MLST + VirulenceFinder: route Listeria samples ---
    lis_inputs = assemblies.contigs.join(kraken_out.kraken_reports)
                  .filter { id, c, rep -> decideOrganism(rep.toString()) == 'listeria' }
                  .map    { id, c, rep -> tuple(id, c) }

    lis_inputs
        .map { id, c -> "MLST + VIRULENCEFINDER will run on: ${id}" }
        .view()

    MLST(lis_inputs)
    VIRULENCEFINDER(lis_inputs)

    // --- Virus / other: no typing, just log it ---
    assemblies.contigs.join(kraken_out.kraken_reports)
        .filter { id, c, rep -> decideOrganism(rep.toString()) == 'virus' }
        .map { id, c, rep -> "VIRUS detected (no typing): ${id}" }
        .view()
}
