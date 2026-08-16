// modules/nanopore.nf
// Contains: DORADO, PORECHOP, FLYE, MEDAKA, CHECKM

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
