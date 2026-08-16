// modules/hybrid.nf
// Contains: UNICYCLER

process UNICYCLER {
    tag "$sample_id"
    publishDir "${params.outdir}/hybrid_assembly", mode: 'copy'
    errorStrategy 'ignore'
    maxRetries 0
    cpus   { params.cpus_spades ?: params.global_cpus }
    memory { params.mem_spades  ?: params.global_memory }

    input:
    tuple val(sample_id), path(illumina_R1), path(illumina_R2), path(ont_reads)

    output:
    tuple val(sample_id), path("${sample_id}_hybrid_contigs.fasta"), optional: true, emit: hybrid_contigs

    script:
    """
    #!/bin/bash
    # Decompress ONT reads if gzipped
    if [[ "${ont_reads}" == *.gz ]]; then
        gunzip -c ${ont_reads} > ont_input.fastq
        ONT_FILE="ont_input.fastq"
    else
        ONT_FILE="${ont_reads}"
    fi

    # Run Unicycler hybrid assembly
    unicycler \\
        -1 ${illumina_R1} \\
        -2 ${illumina_R2} \\
        -l \$ONT_FILE \\
        -o ${sample_id}_unicycler_out \\
        -t ${task.cpus} 2>&1 | tee unicycler_run.log || true

    if [ -f ${sample_id}_unicycler_out/assembly.fasta ]; then
        cp ${sample_id}_unicycler_out/assembly.fasta ${sample_id}_hybrid_contigs.fasta
    else
        echo "Unicycler failed for ${sample_id}" > ${sample_id}_unicycler_failed.txt
        touch ${sample_id}_hybrid_contigs.fasta
    fi
    """
}
