// modules/assembly.nf
// Contains: SPADES, FILTER_CONTIGS, QUAST, BUSCO, BUSCO_SUMMARY

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
