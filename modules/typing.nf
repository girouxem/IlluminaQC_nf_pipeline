// modules/typing.nf
// Contains: KRAKEN2, TB_PROFILER, TB_SUMMARY, SISTR2, SISTR_SUMMARY, AMRFINDER, AMR_SUMMARY, MLST, MLST_SUMMARY, VIRULENCEFINDER, VIRULENCE_SUMMARY, VSNP, VSNP_SUMMARY

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
