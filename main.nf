#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.assembly_mode  = params.assembly_mode  ?: 'isolate'
params.run_nanopore   = params.run_nanopore   ?: false
params.nanopore_input = params.nanopore_input ?: "data_nanopore/*"
params.reads          = params.reads          ?: "data/*_{r1,R1,r2,R2}.fastq.gz"
params.outdir         = params.outdir         ?: "${baseDir}/results"
params.busco_lineage  = params.busco_lineage  ?: 'bacteria_odb10'
params.min_contig_len = params.min_contig_len ?: 500
params.kraken_db      = params.kraken_db      ?: null
params.virulence_db   = params.virulence_db   ?: null
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
params.vsnp_ref       = params.vsnp_ref       ?: null
params.vsnp_options   = params.vsnp_options   ?: null
params.tbprofiler_db  = params.tbprofiler_db  ?: "/home/girouxeml/databases/tbprofiler/tbdb"


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
    path "retention_${sample_id}.csv", emit: retention_csv
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


// ---------- Retention Visualization (always emits files) ----------
process RETENTION_VIZ {
    publishDir "${params.outdir}/retention_plots", mode: 'copy'
    cpus   { params.global_cpus }
    memory { params.global_memory }

    input:
    path retention_csvs

    output:
    path "retention_*.*", emit: retention_files

    script:
    """
    #!/usr/bin/env python3
    import csv, os, sys
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import pandas as pd

    rows = []
    for f in "${retention_csvs}".split():
        if os.path.exists(f):
            df = pd.read_csv(f)
            rows.append(df)
    if not rows:
        with open("retention_summary.tsv", "w") as out:
            out.write("sample\\treads_before\\treads_after\\tbases_before\\tbases_after\\tread_retention_pct\\tbase_retention_pct\\n")
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No retention data", ha='center', va='center')
        plt.savefig("retention_read_bar.png", dpi=150)
        plt.savefig("retention_base_bar.png", dpi=150)
        sys.exit(0)

    data = pd.concat(rows)
    data['read_retention_pct'] = (data['reads_after'] / data['reads_before']) * 100
    data['base_retention_pct'] = (data['bases_after'] / data['bases_before']) * 100
    data.to_csv("retention_summary.tsv", sep="\\t", index=False)

    fig, ax = plt.subplots(figsize=(10, 6))
    x = range(len(data))
    ax.bar(x, data['read_retention_pct'], color='steelblue', label='Read Retention %')
    ax.axhline(y=85, color='red', linestyle='--', label='QC Checkpoint (85%)')
    ax.set_xticks(x)
    ax.set_xticklabels(data['sample'], rotation=45, ha='right')
    ax.set_ylabel('Retention (%)')
    ax.set_title('Read Retention: Raw vs Trimmed')
    ax.legend()
    ax.set_ylim(0, 105)
    plt.tight_layout()
    plt.savefig("retention_read_bar.png", dpi=150)
    plt.close()

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(x, data['base_retention_pct'], color='darkorange', label='Base Retention %')
    ax.axhline(y=85, color='red', linestyle='--', label='QC Checkpoint (85%)')
    ax.set_xticks(x)
    ax.set_xticklabels(data['sample'], rotation=45, ha='right')
    ax.set_ylabel('Retention (%)')
    ax.set_title('Base Retention: Raw vs Trimmed')
    ax.legend()
    ax.set_ylim(0, 105)
    plt.tight_layout()
    plt.savefig("retention_base_bar.png", dpi=150)
    plt.close()
    """
}


// ---------- SPAdes Assembly ----------
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

    # Determine assembly flag based on pipeline parameter
    if [ "${params.assembly_mode}" == "meta" ]; then
        ASM_FLAG="--meta"
    else
        ASM_FLAG="--isolate"
    fi

    spades.py \\
        -1 "${R1}" \\
        -2 "${R2}" \\
        \$ASM_FLAG \\
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


// ---------- QUAST Assembly QC ----------
// ---------- QUAST Assembly QC ----------
process QUAST {
    tag "$sample_id"
    publishDir "${params.outdir}/quast", mode: 'copy'
    errorStrategy 'ignore'
    cpus   { params.cpus_quast ?: params.global_cpus }
    memory { params.mem_quast  ?: params.global_memory }
    input:
    tuple val(sample_id), path(contigs)
    output:
    path "${sample_id}_quast", optional: true, emit: quast_reports
    script:
    """
    #!/bin/bash
    # If assembly is empty or missing, skip QUAST
    if [[ ! -s "${contigs}" ]]; then
        mkdir -p ${sample_id}_quast
        echo "Empty assembly for ${sample_id}" > ${sample_id}_quast/report.txt
        exit 0
    fi

    quast.py ${contigs} -o ${sample_id}_quast -t ${task.cpus} --min-contig 100
    """
}



// ---------- BUSCO Completeness ----------
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
    path "busco_summary.tsv", emit: busco_summary
    script:
    """
    #!/bin/bash
    printf "sample\\tcomplete_pct\\tsingle_pct\\tduplicated_pct\\tfragmented_pct\\tmissing_pct\\ttotal_genes\\n" > busco_summary.tsv
    for d in ${busco_dirs}; do
        sample=\$(basename \$d | sed 's/_busco//')
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
    kraken2 --db ${params.kraken_db} \\
            --memory-mapping \\
            --threads ${task.cpus} \\
            --report ${sample_id}_kraken2.report \\
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

    tb-profiler profile \\
        --read1 ${R1} \\
        --read2 ${R2} \\
        --prefix ${sample_id} \\
        --db ${params.tbprofiler_db} \\
        --txt \\
        --threads ${task.cpus} 2>&1 | tee tbprofiler_run_${sample_id}.log || true

    LINEAGE="NA"
    SUBLINEAGE="NA"
    DRUGS="None"

    if [[ -f results/${sample_id}.results.json ]]; then
        cp results/${sample_id}.results.json . 2>/dev/null || true
        cp results/${sample_id}.results.txt . 2>/dev/null || true

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
    input:
    path tb_tsvs
    output:
    path "tbprofiler_summary.tsv", emit: tb_summary
    script:
    """
    #!/bin/bash
    printf "Sample\\tLineage\\tSublineage\\tDrug_Resistance\\n" > tbprofiler_summary.tsv
    for f in ${tb_tsvs}; do
        if [ -f "\$f" ]; then
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
    sistr -i ${contigs} ${sample_id} -o ${sample_id}_sistr2 -f json 2>&1 | tee sistr_run_${sample_id}.log || true
    if [ -f ${sample_id}_sistr2_predictions.json ]; then
        cp ${sample_id}_sistr2_predictions.json ${sample_id}_sistr2.json
    fi
    if [ -f ${sample_id}_sistr2.json ]; then
        python3 -c "
import json
data=json.load(open('${sample_id}_sistr2.json'))
d = data[0] if isinstance(data, list) and data else {}
serovar = d.get('serovar', d.get('serovar_predicted', d.get('serovar_antigenic', 'NA')))
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
    input:
    path sistr_tsvs
    output:
    path "sistr2_summary.tsv", emit: sistr_summary
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
    input:
    path amr_tsvs
    output:
    path "amrfinder_summary.tsv", emit: amr_summary
    script:
    """
    #!/bin/bash
    printf "Sample\\tAMR_Gene_Count\\tDrug_Classes\\tGene_Symbols\\n" > amrfinder_summary.tsv
    for f in ${amr_tsvs}; do
        if [ -f "\$f" ]; then
            f_name=\$(basename "\$f")
            sample="\${f_name%_amrfinder.tsv}"
            count=\$(tail -n +2 "\$f" | wc -l)
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


// ---------- Combine MLST results into one table ----------
process MLST_SUMMARY {
    publishDir "${params.outdir}/mlst", mode: 'copy'
    input:
    path mlst_tsvs
    output:
    path "mlst_typing_summary.tsv", emit: mlst_summary
    script:
    """
    #!/bin/bash
    printf "Sample\\tMLST_Scheme\\tST_Type\\tAlleles\\n" > mlst_typing_summary.tsv
    for f in ${mlst_tsvs}; do
        if [ -f "\$f" ]; then
            f_name=\$(basename "\$f")
            sample="\${f_name%_mlst.tsv}"
            # mlst output: filepath scheme ST allele1 allele2 ...
            # Extract scheme (col 2), ST (col 3), and join alleles (col 4+) with ;
            scheme=\$(cut -f2 "\$f" | head -1)
            st=\$(cut -f3 "\$f" | head -1)
            alleles=\$(cut -f4- "\$f" | head -1 | tr '\\t' ';')
            printf "%s\\t%s\\t%s\\t%s\\n" "\$sample" "\${scheme:-NA}" "\${st:-NA}" "\${alleles:-NA}" >> mlst_typing_summary.tsv
        fi
    done
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

    # genomicepidemiology/virulencefinder uses -ifa, -d, -o
    virulencefinder.py \\
        -ifa ${contigs} \\
        -o ${sample_id}_vf \\
        -d ${params.virulence_db} \\
        -x \\
        -j 2>&1 | tee virulencefinder_run_${sample_id}.log || true

    # Find and copy the JSON output
    find ${sample_id}_vf -name "*.json" -exec cp {} ${sample_id}_virulence.json \\; 2>/dev/null || echo "{}" > ${sample_id}_virulence.json

    # Find and copy the text output
    find ${sample_id}_vf -name "*.txt" -exec cp {} ${sample_id}_virulence.txt \\; 2>/dev/null || echo "No virulence genes found" > ${sample_id}_virulence.txt
    """
}


// ---------- Combine VirulenceFinder results into one table ----------
process VIRULENCE_SUMMARY {
    publishDir "${params.outdir}/virulencefinder", mode: 'copy'
    input:
    path vir_jsons
    output:
    path "virulencefinder_summary.tsv", emit: vir_summary
    script:
    """
    #!/bin/bash
    printf "Sample\\tVirulence_Gene_Count\\tVirulence_Genes\\n" > virulencefinder_summary.tsv
    for f in ${vir_jsons}; do
        if [ -f "\$f" ] && [ -s "\$f" ]; then
            f_name=\$(basename "\$f")
            sample="\${f_name%_virulence.json}"
            python3 -c "
import json, sys
try:
    d=json.load(open('\$f'))
    genes=[]
    # Try different JSON structures
    if isinstance(d, dict):
        results = d.get('virulencefinder', {}).get('results', {})
        if not results:
            results = d.get('results', {})
        for k,v in results.items():
            if isinstance(v, dict):
                for hit in v.get('virulencefinder', {}).get('hits', []):
                    genes.append(hit.get('gene', 'NA'))
                for hit in v.get('hits', []):
                    genes.append(hit.get('gene', 'NA'))
    count = len(genes)
    genes_str = ';'.join(sorted(set(genes))) if genes else 'None'
    print(f'{count}')
    print(f'{genes_str}')
except Exception as e:
    print('0')
    print('None')
" > /tmp/vf_out.txt 2>/dev/null
            cnt=\$(head -1 /tmp/vf_out.txt)
            genes=\$(tail -1 /tmp/vf_out.txt)
            printf "%s\\t%s\\t%s\\n" "\$sample" "\${cnt:-0}" "\${genes:-None}" >> virulencefinder_summary.tsv
        fi
    done
    """
}

// ---------- MULTIQC ----------
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

// ---------- Checkpoint Report (95/85/80% thresholds) ----------
process CHECKPOINT_REPORT {
    publishDir "${params.outdir}/checkpoints", mode: 'copy'
    cpus   { params.global_cpus }
    memory { params.global_memory }

    input:
    path retention_csvs
    path quast_dirs

    output:
    path "checkpoint_summary.tsv", emit: checkpoint_tsv
    path "checkpoint_bar.png", emit: checkpoint_plot

    script:
    """
    #!/usr/bin/env python3
    import csv, os, sys, json
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import pandas as pd

    # Parse retention CSVs
    retention = {}
    for f in "${retention_csvs}".split():
        if os.path.exists(f):
            df = pd.read_csv(f)
            for _, row in df.iterrows():
                sample = row['sample']
                retention[sample] = {
                    'reads_before': int(row['reads_before']),
                    'reads_after': int(row['reads_after']),
                    'bases_before': int(row['bases_before']),
                    'bases_after': int(row['bases_after'])
                }

    # Parse QUAST reports for assembly stats
    assembly = {}
    for d in "${quast_dirs}".split():
        if os.path.isdir(d):
            report_file = os.path.join(d, 'report.tsv')
            if not os.path.exists(report_file):
                report_file = os.path.join(d, 'report.txt')
            if os.path.exists(report_file):
                with open(report_file) as rf:
                    for line in rf:
                        parts = line.strip().split('\\t')
                        if len(parts) >= 2:
                            if parts[0] == 'Total length':
                                # QUAST dir name is like MBWGS001_quast
                                sample = os.path.basename(d).replace('_quast', '')
                                assembly[sample] = int(parts[1].replace(',', '').replace('.', ''))

    # Build checkpoint summary
    rows = []
    for sample in sorted(set(list(retention.keys()) + list(assembly.keys()))):
        ret = retention.get(sample, {})
        asm = assembly.get(sample, 0)

        reads_before = ret.get('reads_before', 0)
        reads_after = ret.get('reads_after', 0)
        bases_before = ret.get('bases_before', 0)
        bases_after = ret.get('bases_after', 0)

        # Checkpoint 1: Basecalled (always 100% since we start from FASTQs)
        basecalled_pct = 100.0
        basecalled_pass = "PASS" if basecalled_pct >= 95 else "FAIL"

        # Checkpoint 2: QC retention
        qc_pct = (reads_after / reads_before * 100) if reads_before > 0 else 0
        qc_pass = "PASS" if qc_pct >= 85 else "FAIL"

        # Checkpoint 3: Assembly retention (assembly bases / trimmed bases)
        assembly_pct = (asm / bases_after * 100) if bases_after > 0 and asm > 0 else 0
        assembly_pass = "PASS" if assembly_pct >= 80 else "FAIL"

        overall = "PASS" if all([basecalled_pass == "PASS", qc_pass == "PASS", assembly_pass == "PASS"]) else "FAIL"

        rows.append({
            'Sample': sample,
            'Basecalled_Pct': f"{basecalled_pct:.1f}",
            'Basecalled_Status': basecalled_pass,
            'QC_Pct': f"{qc_pct:.1f}",
            'QC_Status': qc_pass,
            'Assembly_Pct': f"{assembly_pct:.1f}" if assembly_pct > 0 else "N/A",
            'Assembly_Status': assembly_pass,
            'Overall_Status': overall
        })

    # Write TSV
    df_out = pd.DataFrame(rows)
    df_out.to_csv("checkpoint_summary.tsv", sep="\\t", index=False)

    # Create bar chart
    if rows:
        fig, ax = plt.subplots(figsize=(12, 6))
        samples = [r['Sample'] for r in rows]
        x = range(len(samples))
        width = 0.25

        basecalled_vals = [float(r['Basecalled_Pct']) for r in rows]
        qc_vals = [float(r['QC_Pct']) for r in rows]
        assembly_vals = [float(r['Assembly_Pct']) if r['Assembly_Pct'] != 'N/A' else 0 for r in rows]

        ax.bar([i - width for i in x], basecalled_vals, width, color='green', label='Basecalled (95%)')
        ax.bar(x, qc_vals, width, color='steelblue', label='QC (85%)')
        ax.bar([i + width for i in x], assembly_vals, width, color='darkorange', label='Assembly (80%)')

        ax.axhline(y=95, color='green', linestyle='--', alpha=0.5)
        ax.axhline(y=85, color='steelblue', linestyle='--', alpha=0.5)
        ax.axhline(y=80, color='darkorange', linestyle='--', alpha=0.5)

        ax.set_xticks(x)
        ax.set_xticklabels(samples, rotation=45, ha='right')
        ax.set_ylabel('Retention (%)')
        ax.set_title('Sample Checkpoints: 95% Basecalled / 85% QC / 80% Assembly')
        ax.legend()
        ax.set_ylim(0, 110)
        plt.tight_layout()
        plt.savefig("checkpoint_bar.png", dpi=150)
        plt.close()
    """
}

// ---------- Assembly Visualization (N50 + Sankey) ----------
process ASSEMBLY_VIZ {
    publishDir "${params.outdir}/assembly_plots", mode: 'copy'
    cpus   { params.global_cpus }
    memory { params.global_memory }

    input:
    path retention_csvs
    path quast_dirs

    output:
    path "n50_bar.png", emit: n50_plot
    path "sankey_retention.png", emit: sankey_plot
    path "assembly_metrics_summary.tsv", emit: assembly_summary

    script:
    """
    #!/usr/bin/env python3
    import csv, os, sys
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import pandas as pd
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    # Parse retention CSVs
    retention = {}
    for f in "${retention_csvs}".split():
        if os.path.exists(f):
            df = pd.read_csv(f)
            for _, row in df.iterrows():
                retention[row['sample']] = {
                    'reads_before': int(row['reads_before']),
                    'reads_after': int(row['reads_after']),
                    'bases_before': int(row['bases_before']),
                    'bases_after': int(row['bases_after'])
                }

    # Parse QUAST reports for N50 and total length
    assembly = {}
    for d in "${quast_dirs}".split():
        if os.path.isdir(d):
            report_file = os.path.join(d, 'report.tsv')
            if not os.path.exists(report_file):
                report_file = os.path.join(d, 'report.txt')
            if os.path.exists(report_file):
                sample = os.path.basename(d).replace('_quast', '')
                with open(report_file) as rf:
                    for line in rf:
                        parts = line.strip().split('\\t')
                        if len(parts) >= 2:
                            if parts[0] == 'N50':
                                assembly[sample] = {'n50': int(parts[1].replace(',', '').replace('.', ''))}
                            if parts[0] == 'Total length' and sample in assembly:
                                assembly[sample]['total'] = int(parts[1].replace(',', '').replace('.', ''))

    # ---- N50 Bar Chart ----
    samples = sorted(assembly.keys())
    if samples:
        fig, ax = plt.subplots(figsize=(12, 6))
        n50_vals = [assembly[s].get('n50', 0) for s in samples]
        colors = ['steelblue' if v > 0 else 'red' for v in n50_vals]
        ax.bar(range(len(samples)), n50_vals, color=colors)
        ax.set_xticks(range(len(samples)))
        ax.set_xticklabels(samples, rotation=45, ha='right')
        ax.set_ylabel('N50 (bp)')
        ax.set_title('Assembly N50 per Sample')
        plt.tight_layout()
        plt.savefig("n50_bar.png", dpi=150)
        plt.close()
    else:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No assembly data", ha='center')
        plt.savefig("n50_bar.png", dpi=150)
        plt.close()

    # ---- Summary TSV ----
    rows = []
    for s in sorted(set(list(retention.keys()) + list(assembly.keys()))):
        ret = retention.get(s, {})
        asm = assembly.get(s, {})
        rows.append({
            'Sample': s,
            'Reads_Before': ret.get('reads_before', 0),
            'Reads_After': ret.get('reads_after', 0),
            'Bases_Before': ret.get('bases_before', 0),
            'Bases_After': ret.get('bases_after', 0),
            'Assembly_Length': asm.get('total', 0),
            'N50': asm.get('n50', 0),
            'Read_Retention_Pct': f"{(ret.get('reads_after',0)/ret.get('reads_before',1)*100):.1f}" if ret.get('reads_before',0) > 0 else "0",
            'Assembly_Retention_Pct': f"{(asm.get('total',0)/ret.get('bases_after',1)*100):.1f}" if ret.get('bases_after',0) > 0 and asm.get('total',0) > 0 else "0"
        })
    pd.DataFrame(rows).to_csv("assembly_metrics_summary.tsv", sep="\\t", index=False)

    # ---- Sankey Diagram ----
    # ---- Retention Flow Chart (Matplotlib stacked bar) ----
    if rows:
        all_samples = [r['Sample'] for r in rows]
        reads_before = [r['Reads_Before'] for r in rows]
        reads_after = [r['Reads_After'] for r in rows]
        assembly_reads = [r['Assembly_Length'] for r in rows]

        fig, ax = plt.subplots(figsize=(12, 6))
        x = range(len(all_samples))
        width = 0.25

        ax.bar([i - width for i in x], reads_before, width, color='steelblue', label='Raw Reads')
        ax.bar(x, reads_after, width, color='darkorange', label='Trimmed Reads')
        ax.bar([i + width for i in x], assembly_reads, width, color='green', label='Assembly Bases')

        ax.set_xticks(x)
        ax.set_xticklabels(all_samples, rotation=45, ha='right')
        ax.set_ylabel('Count (Reads / Bases)')
        ax.set_title('Retention Flow: Raw → Trimmed → Assembly')
        ax.legend()
        ax.set_yscale('log') # Log scale helps visualize vastly different magnitudes
        plt.tight_layout()
        plt.savefig("sankey_retention.png", dpi=150)
        plt.close()
    else:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No data for flow chart", ha='center')
        plt.savefig("sankey_retention.png", dpi=150)
        plt.close()

    """
}
// ---------- Read Length Histograms ----------
process READ_LEN_HIST {
    publishDir "${params.outdir}/read_length_plots", mode: 'copy'
    cpus   { params.global_cpus }
    memory { params.global_memory }

    input:
    tuple val(sample_id), path(raw_R1), path(raw_R2), path(trim_R1), path(trim_R2)

    output:
    tuple val(sample_id), path("readlen_${sample_id}.png"), emit: hist_plot

    script:
    """
    #!/usr/bin/env python3
    import gzip, os, sys
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    def get_read_lengths(fastq):
        lengths = []
        if fastq.endswith('.gz'):
            fh = gzip.open(fastq, 'rt')
        else:
            fh = open(fastq, 'r')
        count = 0
        for i, line in enumerate(fh):
            if i % 4 == 1:  # sequence line
                lengths.append(len(line.strip()))
                count += 1
                if count >= 10000:  # sample first 10k reads for speed
                    break
        fh.close()
        return lengths

    raw_lens = get_read_lengths("${raw_R1}")
    trim_lens = get_read_lengths("${trim_R1}")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(raw_lens, bins=50, alpha=0.6, color='steelblue', label='Raw')
    ax.hist(trim_lens, bins=50, alpha=0.6, color='darkorange', label='Trimmed')
    ax.set_xlabel('Read Length (bp)')
    ax.set_ylabel('Count')
    ax.set_title('Read Length Distribution: ${sample_id}')
    ax.legend()
    plt.tight_layout()
    plt.savefig("readlen_${sample_id}.png", dpi=150)
    plt.close()
    """
}

// ---------- vSNP3 (M. bovis phylogenetic grouping) ----------
process VSNP {
    tag "$sample_id"
    publishDir "${params.outdir}/vsnp", mode: 'copy'
    errorStrategy 'ignore'
    maxRetries 0
    cpus   { params.global_cpus }
    memory { params.global_memory }

    input:
    tuple val(sample_id), path(R1), path(R2)

    output:
    tuple val(sample_id), path("${sample_id}_vsnp_group.tsv"), optional: true, emit: vsnp_result
    path "${sample_id}_vsnp_failed.txt", optional: true

    when:
    params.vsnp_ref != null

    script:
    """
    #!/bin/bash

    if [[ ! -s "${R1}" || ! -s "${R2}" ]]; then
        echo "No reads for ${sample_id}" > ${sample_id}_vsnp_failed.txt
        printf "Sample\\tvSNP_Group\\n" > ${sample_id}_vsnp_group.tsv
        printf "${sample_id}\\tNA\\n" >> ${sample_id}_vsnp_group.tsv
        exit 0
    fi

    # ============ STEP 1: Align reads and call SNPs ============
    vsnp3_step1.py \\
        -r1 ${R1} \\
        -r2 ${R2} \\
        -r ${params.vsnp_ref} \\
        -n ${sample_id} \\
        -o ${sample_id}_vsnp_out 2>&1 | tee vsnp_run_${sample_id}.log || true

    # Find the VCF file produced by step1
    VCF_FILE=\$(find ${sample_id}_vsnp_out -name "*_zc.vcf" -type f 2>/dev/null | head -1)

    if [ -z "\$VCF_FILE" ]; then
        echo "No VCF produced by step1" > ${sample_id}_vsnp_failed.txt
        printf "Sample\\tvSNP_Group\\n" > ${sample_id}_vsnp_group.tsv
        printf "${sample_id}\\tNA\\n" >> ${sample_id}_vsnp_group.tsv
        exit 0
    fi

    echo "=== Found VCF: \$VCF_FILE ===" >> vsnp_run_${sample_id}.log

    # ============ STEP 2: Group assignment using defining SNPs ============
    # Copy VCF to a clean directory for step2
    mkdir -p ${sample_id}_step2_input
    cp "\$VCF_FILE" ${sample_id}_step2_input/

    STEP2_ARGS="-wd ${sample_id}_step2_input -o ${sample_id}_step2_out --show_groups"

    # Add defining SNPs file if provided
    if [ -n "${params.vsnp_options}" ] && [ -f "${params.vsnp_options}" ]; then
        STEP2_ARGS="\$STEP2_ARGS -s ${params.vsnp_options}"
    fi

    echo "=== Running: vsnp3_step2.py \$STEP2_ARGS ===" >> vsnp_run_${sample_id}.log
    vsnp3_step2.py \$STEP2_ARGS 2>&1 | tee vsnp_step2_run_${sample_id}.log || true

        # ============ Parse group assignment from step2 output ============
    # The group name IS the subdirectory name inside step2_out
    group="NA"
    if [ -d ${sample_id}_step2_out ]; then
        echo "=== Step2 output files ===" >> vsnp_run_${sample_id}.log
        find ${sample_id}_step2_out -type f >> vsnp_run_${sample_id}.log 2>/dev/null

        # Extract group names from subdirectories (ignore files and temp folders)
        groups=\$(find ${sample_id}_step2_out -mindepth 1 -maxdepth 1 -type d -printf "%f\\n" 2>/dev/null | grep -v "step2_is_running" | sort -r | tr '\\n' ';')
        # Remove trailing semicolon
        groups="\${groups%;}"
        
        if [ -n "\$groups" ]; then
            group="\$groups"
        fi
    else
        echo "Step2 produced no output directory" >> vsnp_run_${sample_id}.log
    fi

    echo "=== Final group: \$group ===" >> vsnp_run_${sample_id}.log
    printf "Sample\\tvSNP_Group\\n" > ${sample_id}_vsnp_group.tsv
    printf "${sample_id}\\t\${group:-NA}\\n" >> ${sample_id}_vsnp_group.tsv
    """
}

// ---------- Combine vSNP results into one table ----------
process VSNP_SUMMARY {
    publishDir "${params.outdir}/vsnp", mode: 'copy'
    input:
    path vsnp_tsvs
    output:
    path "vsnp_summary.tsv", emit: vsnp_summary
    script:
    """
    #!/bin/bash
    printf "Sample\\tvSNP_Group\\n" > vsnp_summary.tsv
    for f in ${vsnp_tsvs}; do
        if [ -f "\$f" ]; then
            tail -n +2 "\$f" >> vsnp_summary.tsv
        fi
    done
    """
}
// ===================== NANOPORE BRANCH =====================

// ---------- Dorado Basecaller ----------
process DORADO {
    tag "$sample_id"
    publishDir "${params.outdir}/nanopore_basecalled", mode: 'copy'
    errorStrategy 'ignore'
    maxRetries 0
    cpus   { params.global_cpus }
    memory { params.global_memory }

    input:
    tuple val(sample_id), path(raw_signal)

    output:
    tuple val(sample_id), path("${sample_id}_basecalled.fastq.gz"), optional: true, emit: basecalled_fastq

    script:
    """
    #!/bin/bash
    # Dorado basecaller - supports both POD5 and FAST5
    # Auto-detect format by extension
    INPUT_DIR="${sample_id}_input"
    mkdir -p \$INPUT_DIR
    cp ${raw_signal} \$INPUT_DIR/

    # Run basecaller with duplex mode and 5mC modification
    dorado basecaller sup,duplex \$INPUT_DIR > ${sample_id}_basecalled.fastq 2> dorado_run.log || true
    
    # Also try demultiplexing if barcoded
    if [ -s ${sample_id}_basecalled.fastq ]; then
        gzip ${sample_id}_basecalled.fastq
    else
        echo "Dorado basecalling failed for ${sample_id}" > ${sample_id}_dorad_failed.txt
        touch ${sample_id}_basecalled.fastq.gz
    fi
    """
}

// ---------- PoreChop (adapter trimming for ONT) ----------
process PORECHOP {
    tag "$sample_id"
    publishDir "${params.outdir}/nanopore_trimmed", mode: 'copy'
    errorStrategy 'ignore'
    maxRetries 0
    cpus   { params.global_cpus }
    memory '8 GB'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_porechop.fastq.gz"), optional: true, emit: trimmed_fastq

    script:
    """
    #!/bin/bash
    # Decompress input if gzipped, PoreChop handles plain FASTQ better
    if [[ "${reads}" == *.gz ]]; then
        gunzip -c ${reads} > input.fastq
        INPUT_FILE="input.fastq"
    else
        INPUT_FILE="${reads}"
    fi

    # Run PoreChop  no --format flag, let it auto-detect
    porechop -i \$INPUT_FILE -o ${sample_id}_porechop.fastq --threads ${task.cpus} 2>&1 | tee porechop_run.log

    if [ -s ${sample_id}_porechop.fastq ]; then
        gzip ${sample_id}_porechop.fastq
    else
        echo "PoreChop failed for ${sample_id}" > ${sample_id}_porechop_failed.txt
        touch ${sample_id}_porechop.fastq.gz
    fi
    """
}


// ---------- Flye (long-read assembly) ----------
process FLYE {
    tag "$sample_id"
    publishDir "${params.outdir}/nanopore_assembly", mode: 'copy'
    errorStrategy 'ignore'
    maxRetries 0
    cpus   { params.cpus_spades ?: params.global_cpus }
    memory { params.mem_spades  ?: params.global_memory }

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_flye_contigs.fasta"), optional: true, emit: contigs

    script:
    """
    #!/bin/bash
    # Always decompress to plain FASTQ
    if [[ "${reads}" == *.gz ]]; then
        gunzip -c ${reads} > ${sample_id}_input.fastq
        INPUT_FILE="${sample_id}_input.fastq"
    else
        cp ${reads} ${sample_id}_input.fastq
        INPUT_FILE="${sample_id}_input.fastq"
    fi

    # Check if file has content
    if [ ! -s \$INPUT_FILE ]; then
        echo "Empty input file for ${sample_id}" > ${sample_id}_flye_failed.txt
        touch ${sample_id}_flye_contigs.fasta
        exit 0
    fi

    # Determine assembly flag
    if [ "${params.assembly_mode}" == "meta" ]; then
        ASM_FLAG="--meta"
    else
        ASM_FLAG=""
    fi

    # Run Flye  use --nano-raw instead of --nano-hq for noisy metagenomic data
    # Lower minimum read length to 200 (default is 1000)
    flye --nano-raw \$INPUT_FILE \\
         --out-dir ${sample_id}_flye_out \\
         --threads ${task.cpus} \\
         \$ASM_FLAG 2>&1 | tee flye_run.log || true

    if [ -f ${sample_id}_flye_out/assembly.fasta ]; then
        cp ${sample_id}_flye_out/assembly.fasta ${sample_id}_flye_contigs.fasta
    else
        echo "Flye assembly failed for ${sample_id}" > ${sample_id}_flye_failed.txt
        touch ${sample_id}_flye_contigs.fasta
    fi
    """
}


// ---------- Medaka (polishing) ----------
process MEDAKA {
    tag "$sample_id"
    publishDir "${params.outdir}/nanopore_polished", mode: 'copy'
    errorStrategy 'ignore'
    maxRetries 0
    cpus   { params.global_cpus }
    memory { params.global_memory }

    input:
    tuple val(sample_id), path(assembly), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_polished.fasta"), optional: true, emit: polished_contigs

    script:
    """
    #!/bin/bash
    # Decompress reads if gzipped
    if [[ "${reads}" == *.gz ]]; then
        gunzip -c ${reads} > ${sample_id}_reads.fastq
        READS="${sample_id}_reads.fastq"
    else
        READS="${reads}"
    fi

    medaka_consensus -i \$READS \\
                     -d ${assembly} \\
                     -o ${sample_id}_medaka_out \\
                     -t ${task.cpus} 2>&1 | tee medaka_run.log || true

    if [ -f ${sample_id}_medaka_out/consensus.fasta ]; then
        cp ${sample_id}_medaka_out/consensus.fasta ${sample_id}_polished.fasta
    else
        # If medaka fails, use the unpolished assembly
        cp ${assembly} ${sample_id}_polished.fasta
        echo "Medaka polishing failed for ${sample_id}  using unpolished assembly" > ${sample_id}_medaka_failed.txt
    fi
    """
}

// ---------- CheckM (bacterial completeness + contamination) ----------
// ---------- CheckM (bacterial completeness + contamination) ----------
process CHECKM {
    tag "$sample_id"
    publishDir "${params.outdir}/nanopore_checkm", mode: 'copy'
    errorStrategy 'ignore'
    maxRetries 0
    cpus   { params.global_cpus }
    memory { params.global_memory }

    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), path("${sample_id}_checkm.tsv"), optional: true, emit: checkm_result

    script:
    """
    #!/bin/bash
    if [[ ! -s "${contigs}" ]]; then
        printf "Bin Id\\tMarker lineage\\t# genomes\\t# markers\\t# marker sets\\t0\\t1\\t2\\t3\\t4\\t5+\\tCompleteness\\tContamination\\tStrain heterogeneity\\n" > ${sample_id}_checkm.tsv
        printf "${sample_id}\\tNA\\t0\\t0\\t0\\t0\\t0\\t0\\t0\\t0\\t0\\tNA\\tNA\\tNA\\n" >> ${sample_id}_checkm.tsv
        exit 0
    fi

    mkdir -p ${sample_id}_checkm_input
    cp ${contigs} ${sample_id}_checkm_input/
    
    checkm lineage_wf ${sample_id}_checkm_input ${sample_id}_checkm_out -t ${task.cpus} --file ${sample_id}_checkm.tsv 2>&1 | tee checkm_run.log || true

    if [ ! -s ${sample_id}_checkm.tsv ]; then
        printf "Bin Id\\tMarker lineage\\t# genomes\\t# markers\\t# marker sets\\t0\\t1\\t2\\t3\\t4\\t5+\\tCompleteness\\tContamination\\tStrain heterogeneity\\n" > ${sample_id}_checkm.tsv
        printf "${sample_id}\\tNA\\t0\\t0\\t0\\t0\\t0\\t0\\t0\\t0\\t0\\tNA\\tNA\\tNA\\n" >> ${sample_id}_checkm.tsv
    fi
    """
}



// Full pipeline
workflow {
    if (params.run_nanopore) {
        nanopore_pipeline()
    } else {
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
        trimmed_fastqc = FASTQC_TRIMMED(trimmed.trimmed_reads)

        retention_input = raw_pairs.join(trimmed.trimmed_reads)
        RETENTION_TRACK(retention_input)
        retention_csvs = RETENTION_TRACK.out.retention_csv.collect()
        RETENTION_VIZ(retention_csvs)

        assemblies = SPADES(trimmed.trimmed_reads)
        filtered   = FILTER_CONTIGS(assemblies.contigs)
        quast_out  = QUAST(assemblies.contigs)
        busco_out  = BUSCO(filtered.filtered_contigs)
        busco_collected = busco_out.busco_reports.collect()
        BUSCO_SUMMARY(busco_collected)

        kraken_out = KRAKEN2(assemblies.contigs)

        vsnp_summary_file = channel.of()

        if (params.kraken_db) {
            tb_inputs = trimmed.trimmed_reads.join(kraken_out.kraken_reports)
                         .filter { id, r1, r2, rep -> decideOrganism(rep.toString()) == 'mbovis' }
                         .map    { id, r1, r2, rep -> tuple(id, r1, r2) }
            TB_PROFILER(tb_inputs)
            tb_mqc_reports = TB_PROFILER.out.tbprofiler_mqc.map { id, f -> f }.collect()
            TB_SUMMARY(tb_mqc_reports)

            if (params.vsnp_ref) {
                VSNP(tb_inputs)
                vsnp_tsv = VSNP.out.vsnp_result.map { id, f -> f }.collect()
                VSNP_SUMMARY(vsnp_tsv)
                vsnp_summary_file = VSNP_SUMMARY.out.vsnp_summary.collect().flatten()
            }

            sal_inputs = assemblies.contigs.join(kraken_out.kraken_reports)
                          .filter { id, c, rep -> decideOrganism(rep.toString()) == 'salmonella' }
                          .map    { id, c, rep -> tuple(id, c) }
            SISTR2(sal_inputs)
            sistr_tsv = SISTR2.out.sistr2_tsv.map { id, f -> f }.collect()
            SISTR_SUMMARY(sistr_tsv)

            AMRFINDER(sal_inputs)
            amr_reports = AMRFINDER.out.amr_tsv.map { id, f -> f }.collect()
            AMR_SUMMARY(amr_reports)

            lis_inputs = assemblies.contigs.join(kraken_out.kraken_reports)
                          .filter { id, c, rep -> decideOrganism(rep.toString()) == 'listeria' }
                          .map    { id, c, rep -> tuple(id, c) }
            MLST(lis_inputs)
            mlst_tsv = MLST.out.mlst_tsv.map { id, f -> f }.collect()
            MLST_SUMMARY(mlst_tsv)

            VIRULENCEFINDER(lis_inputs)
            vir_jsons = VIRULENCEFINDER.out.vir_json.map { id, f -> f }.collect()
            VIRULENCE_SUMMARY(vir_jsons)
        }

        raw_reports     = raw_fastqc.fastqc_zip.map { id, f -> f }
        trimmed_reports = trimmed.fastp_json.map { id, f -> f }
        post_reports    = trimmed_fastqc.fastqc_zip_trim.map { id, f -> f }
        quast_reports   = quast_out.quast_reports.map { p -> p }
        busco_reports   = busco_out.busco_reports.map { p -> p }

        retention_files = RETENTION_VIZ.out.retention_files.collect().flatten()
        tb_summary_file = TB_SUMMARY.out.tb_summary.collect().flatten()
        sistr_summary_file = SISTR_SUMMARY.out.sistr_summary.collect().flatten()
        amr_summary_file = AMR_SUMMARY.out.amr_summary.collect().flatten()
        busco_summary_file = BUSCO_SUMMARY.out.busco_summary.collect().flatten()
        mlst_summary_file = MLST_SUMMARY.out.mlst_summary.collect().flatten()
        vir_summary_file  = VIRULENCE_SUMMARY.out.vir_summary.collect().flatten()

        all_reports = raw_reports
            .mix(trimmed_reports, post_reports, quast_reports, busco_reports)
            .mix(retention_files)
            .mix(tb_summary_file)
            .mix(vsnp_summary_file)
            .mix(sistr_summary_file)
            .mix(amr_summary_file)
            .mix(busco_summary_file)
            .mix(mlst_summary_file)
            .mix(vir_summary_file)
            .collect()

        MULTIQC(all_reports, file("${baseDir}/multiqc_config.yaml"))
    }
}




workflow test_SANITY_CHECK {
    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    SANITY_CHECK(raw_pairs)
}

workflow test_FASTQC_RAW {
    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    SANITY_CHECK(raw_pairs)
    FASTQC_RAW(raw_pairs)
}

workflow test_FASTP {
    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    SANITY_CHECK(raw_pairs)
    FASTP(raw_pairs)
}

workflow test_FASTQC_TRIMMED {
    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    trimmed = FASTP(raw_pairs)
    trimmed_for_qc = trimmed.trimmed_reads.map { id, r1, r2 -> tuple(id, r1, r2) }
    FASTQC_TRIMMED(trimmed_for_qc)
}

workflow test_RETENTION_TRACK {
    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    trimmed = FASTP(raw_pairs)
    retention_input = raw_pairs.join(trimmed.trimmed_reads)
    RETENTION_TRACK(retention_input)
}

workflow test_SPADES {
    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    trimmed = FASTP(raw_pairs)
    SPADES(trimmed.trimmed_reads)
}

workflow test_QUAST {
    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    trimmed = FASTP(raw_pairs)
    assemblies = SPADES(trimmed.trimmed_reads)
    QUAST(assemblies.contigs)
}

workflow test_BUSCO {
    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    trimmed = FASTP(raw_pairs)
    assemblies = SPADES(trimmed.trimmed_reads)
    filtered = FILTER_CONTIGS(assemblies.contigs)
    busco_out = BUSCO(filtered.filtered_contigs)
    busco_collected = busco_out.busco_reports.collect()
    BUSCO_SUMMARY(busco_collected)
}

workflow test_KRAKEN2 {
    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    trimmed = FASTP(raw_pairs)
    assemblies = SPADES(trimmed.trimmed_reads)
    KRAKEN2(assemblies.contigs)
}

workflow test_TB_PROFILER {
    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    trimmed = FASTP(raw_pairs)
    TB_PROFILER(trimmed.trimmed_reads)
}

workflow test_KRAKEN2_ROUTING {
    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    trimmed = FASTP(raw_pairs)
    assemblies = SPADES(trimmed.trimmed_reads)
    kraken_out = KRAKEN2(assemblies.contigs)

    kraken_out.kraken_reports.map { id, rep -> "KRAKEN2: ${id} -> ${rep}" }.view()

    routed = assemblies.contigs.join(kraken_out.kraken_reports)
    routed.map { id, contigs, rep ->
        def org = decideOrganism(rep.toString())
        "ROUTING: ${id} -> ${org}"
    }.view()

    tb_inputs = trimmed.trimmed_reads.join(kraken_out.kraken_reports)
                 .filter { id, r1, r2, rep -> decideOrganism(rep.toString()) == 'mbovis' }
                 .map    { id, r1, r2, rep -> tuple(id, r1, r2) }
    tb_inputs.map { id, r1, r2 -> "TB_PROFILER will run on: ${id}" }.view()
    TB_PROFILER(tb_inputs)

    sal_inputs = assemblies.contigs.join(kraken_out.kraken_reports)
                  .filter { id, c, rep -> decideOrganism(rep.toString()) == 'salmonella' }
                  .map    { id, c, rep -> tuple(id, c) }
    sal_inputs.map { id, c -> "SISTR2 + AMRFINDER will run on: ${id}" }.view()
    SISTR2(sal_inputs)
    AMRFINDER(sal_inputs)

    lis_inputs = assemblies.contigs.join(kraken_out.kraken_reports)
                  .filter { id, c, rep -> decideOrganism(rep.toString()) == 'listeria' }
                  .map    { id, c, rep -> tuple(id, c) }
    lis_inputs.map { id, c -> "MLST + VIRULENCEFINDER will run on: ${id}" }.view()
    MLST(lis_inputs)
    VIRULENCEFINDER(lis_inputs)
}

/*
 * TEST: VirulenceFinder only
 * Run with:
 *   nextflow run . -entry test_VIRULENCEFINDER --virulence_db /home/girouxeml/databases/virulencefinder_db -resume
 */
workflow test_VIRULENCEFINDER {
    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    trimmed = FASTP(raw_pairs)
    assemblies = SPADES(trimmed.trimmed_reads)
    
    assemblies.contigs
        .map { it.toString() }
        .view()
    
    VIRULENCEFINDER(assemblies.contigs)
}

/*
 * TEST: vSNP3 only
 * Run with:
 *   nextflow run . -entry test_VSNP --vsnp_ref /home/girouxeml/databases/vsnp/Mbovis_AF2122.fa -resume
 */
workflow test_VSNP {
    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    trimmed = FASTP(raw_pairs)
    
    trimmed.trimmed_reads
        .map { it.toString() }
        .view()
    
    VSNP(trimmed.trimmed_reads)
}
/*
 * NANOPORE PIPELINE
 * Run with:
 *   nextflow run . --run_nanopore true
 *
 * Input: data_nanopore/ directory with:
 *   - *.pod5 files (raw signal  will be basecalled with Dorado)
 *   - *.fast5 files (raw signal  will be basecalled with Dorado)
 *   - *.fastq.gz files (already basecalled  skips Dorado)
 */
workflow nanopore_pipeline {
    log.info "=== Running Nanopore Pipeline ==="
    log.info "Nanopore input: ${params.nanopore_input}"
    
    // Create input channel from directory
    // Each file or subdirectory becomes one sample
    ont_inputs = channel.fromPath(params.nanopore_input, type: 'file')
        .map { f -> tuple(f.baseName, f) }
    
    // Route based on file extension
    // POD5/FAST5 → Dorado → PoreChop → Flye → Medaka → QUAST + CheckM
    // FASTQ → PoreChop → Flye → Medaka → QUAST + CheckM
    
    needs_basecalling = ont_inputs.filter { id, f -> f.name.endsWith('.pod5') || f.name.endsWith('.fast5') }
    already_basecalled = ont_inputs.filter { id, f -> f.name.endsWith('.fastq') || f.name.endsWith('.fastq.gz') }
    
    // Basecall raw signal
    basecalled = DORADO(needs_basecalling)
    
    // Merge basecalled + already-FASTQ into one channel
    all_ont_reads = basecalled.basecalled_fastq.mix(already_basecalled)
    
    // Trim adapters
    trimmed_ont = PORECHOP(all_ont_reads)
    
    // COMMENTED OUT FOR PORECHOP TESTING
    // Assemble
    ont_assemblies = FLYE(trimmed_ont.trimmed_fastq)
    
    // Polish (needs both assembly + reads)
    polished_inputs = ont_assemblies.contigs.join(trimmed_ont.trimmed_fastq)
    polished_ont = MEDAKA(polished_inputs)
    
    quast_ont_out = QUAST(polished_ont.polished_contigs)
    checkm_ont_out = CHECKM(polished_ont.polished_contigs)
    
    log.info "Nanopore pipeline complete. Check results/nanopore_assembly/ for contigs."
}
