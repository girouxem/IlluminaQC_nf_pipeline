#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ===================== INCLUDE MODULES =====================
include { SANITY_CHECK; FASTQC_RAW; FASTP; FASTQC_TRIMMED } from './modules/qc'
include { RETENTION_TRACK; RETENTION_VIZ } from './modules/retention'
include { SPADES; FILTER_CONTIGS; QUAST; BUSCO; BUSCO_SUMMARY } from './modules/assembly'
include { KRAKEN2; TB_PROFILER; TB_SUMMARY; SISTR2; SISTR_SUMMARY; AMRFINDER; AMR_SUMMARY; MLST; MLST_SUMMARY; VIRULENCEFINDER; VIRULENCE_SUMMARY; VSNP; VSNP_SUMMARY } from './modules/typing'
include { DORADO; PORECHOP; FLYE; MEDAKA; CHECKM } from './modules/nanopore'
include { UNICYCLER } from './modules/hybrid'
include { CHECKPOINT_REPORT; ASSEMBLY_VIZ; READ_LEN_HIST } from './modules/viz'
include { MULTIQC } from './modules/multiqc'

// ===================== PARAMETERS =====================
params.run_hybrid     = params.run_hybrid     ?: false
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

// ===================== HELPER FUNCTIONS =====================
def decideOrganism(String kraken_report_path) {
    def lines = new File(kraken_report_path).readLines()
    def species = lines.find { it.contains("\tS\t") } ?: ''
    if (species.contains("Mycobacterium")) return "mbovis"
    if (species.contains("Salmonella"))    return "salmonella"
    if (species.contains("Listeria"))      return "listeria"
    if (species.toLowerCase().contains("virus")) return "virus"
    return "bacteria"
}

// ===================== MAIN WORKFLOW =====================
workflow {
    if (params.run_nanopore) {
        nanopore_pipeline()
    } else if (params.run_hybrid) {
        hybrid_assembly()
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

// ===================== NANOPORE WORKFLOW =====================
workflow nanopore_pipeline {
    log.info "=== Running Nanopore Pipeline ==="
    log.info "Nanopore input: ${params.nanopore_input}"
    
    ont_inputs = channel.fromPath(params.nanopore_input, type: 'file')
        .map { f -> tuple(f.baseName.replaceFirst(/\.(fastq|fq)(\.gz)?$/, ''), f) }

    needs_basecalling = ont_inputs.filter { id, f -> f.name.endsWith('.pod5') || f.name.endsWith('.fast5') }
    already_basecalled = ont_inputs.filter { id, f -> f.name.endsWith('.fastq') || f.name.endsWith('.fastq.gz') }
    
    basecalled = DORADO(needs_basecalling)
    all_ont_reads = basecalled.basecalled_fastq.mix(already_basecalled)
    
    trimmed_ont = PORECHOP(all_ont_reads)
    ont_assemblies = FLYE(trimmed_ont.trimmed_fastq)
    
    polished_inputs = ont_assemblies.contigs.join(trimmed_ont.trimmed_fastq)
    polished_ont = MEDAKA(polished_inputs)
    
    quast_ont_out = QUAST(polished_ont.polished_contigs)
    checkm_ont_out = CHECKM(polished_ont.polished_contigs)
    
    kraken_ont_out = KRAKEN2(polished_ont.polished_contigs)
    
    tb_summary_ch    = channel.of()
    vsnp_summary_ch  = channel.of()
    sistr_summary_ch = channel.of()
    amr_summary_ch   = channel.of()
    mlst_summary_ch  = channel.of()
    vir_summary_ch   = channel.of()
    
    if (params.kraken_db) {
        tb_ont_inputs = trimmed_ont.trimmed_fastq.join(kraken_ont_out.kraken_reports)
                         .filter { id, fq, rep -> decideOrganism(rep.toString()) == 'mbovis' }
                         .map    { id, fq, rep -> tuple(id, fq, fq) }
        TB_PROFILER(tb_ont_inputs)
        tb_mqc_reports = TB_PROFILER.out.tbprofiler_mqc.map { id, f -> f }.collect()
        TB_SUMMARY(tb_mqc_reports)
        tb_summary_ch = TB_SUMMARY.out.tb_summary
        
        if (params.vsnp_ref) {
            VSNP(tb_ont_inputs)
            vsnp_tsv = VSNP.out.vsnp_result.map { id, f -> f }.collect()
            VSNP_SUMMARY(vsnp_tsv)
            vsnp_summary_ch = VSNP_SUMMARY.out.vsnp_summary.collect().flatten()
        }
        
        sal_ont_inputs = polished_ont.polished_contigs.join(kraken_ont_out.kraken_reports)
                          .filter { id, c, rep -> decideOrganism(rep.toString()) == 'salmonella' }
                          .map    { id, c, rep -> tuple(id, c) }
        SISTR2(sal_ont_inputs)
        sistr_tsv = SISTR2.out.sistr2_tsv.map { id, f -> f }.collect()
        SISTR_SUMMARY(sistr_tsv)
        sistr_summary_ch = SISTR_SUMMARY.out.sistr_summary
        
        AMRFINDER(sal_ont_inputs)
        amr_reports = AMRFINDER.out.amr_tsv.map { id, f -> f }.collect()
        AMR_SUMMARY(amr_reports)
        amr_summary_ch = AMR_SUMMARY.out.amr_summary
        
        lis_ont_inputs = polished_ont.polished_contigs.join(kraken_ont_out.kraken_reports)
                          .filter { id, c, rep -> decideOrganism(rep.toString()) == 'listeria' }
                          .map    { id, c, rep -> tuple(id, c) }
        MLST(lis_ont_inputs)
        mlst_tsv = MLST.out.mlst_tsv.map { id, f -> f }.collect()
        MLST_SUMMARY(mlst_tsv)
        mlst_summary_ch = MLST_SUMMARY.out.mlst_summary
        
        VIRULENCEFINDER(lis_ont_inputs)
        vir_jsons = VIRULENCEFINDER.out.vir_json.map { id, f -> f }.collect()
        VIRULENCE_SUMMARY(vir_jsons)
        vir_summary_ch = VIRULENCE_SUMMARY.out.vir_summary
    }
    
    quast_ont_reports = quast_ont_out.quast_reports.map { p -> p }
    
    all_reports = quast_ont_reports
        .mix(tb_summary_ch, vsnp_summary_ch, sistr_summary_ch, amr_summary_ch, mlst_summary_ch, vir_summary_ch)
        .collect()
    
    MULTIQC(all_reports, file("${baseDir}/multiqc_config.yaml"))
    
    log.info "Nanopore pipeline complete. Check results/nanopore_assembly/ for contigs."
}

// ===================== HYBRID WORKFLOW =====================
workflow hybrid_assembly {
    log.info "=== Running Hybrid Assembly Pipeline ==="
    
    illumina_pairs = channel.fromFilePairs(params.reads, flat: false)
    illumina_trimmed = FASTP(illumina_pairs)
    
    ont_inputs = channel.fromPath(params.nanopore_input, type: 'file')
        .map { f -> tuple(f.baseName.replaceFirst(/\.(fastq|fq)(\.gz)?$/, ''), f) }

    needs_basecalling = ont_inputs.filter { id, f -> f.name.endsWith('.pod5') || f.name.endsWith('.fast5') }
    already_basecalled = ont_inputs.filter { id, f -> f.name.endsWith('.fastq') || f.name.endsWith('.fastq.gz') }
    basecalled = DORADO(needs_basecalling)
    all_ont_reads = basecalled.basecalled_fastq.mix(already_basecalled)
    ont_trimmed = PORECHOP(all_ont_reads)
    
    hybrid_inputs = illumina_trimmed.trimmed_reads.join(ont_trimmed.trimmed_fastq)
    
    hybrid_inputs.map { id, r1, r2, ont -> "HYBRID: ${id} has Illumina + Nanopore" }.view()
    
    hybrid_assemblies = UNICYCLER(hybrid_inputs)
    
    quast_hybrid_out = QUAST(hybrid_assemblies.hybrid_contigs)
    checkm_hybrid_out = CHECKM(hybrid_assemblies.hybrid_contigs)
    
    kraken_hybrid_out = KRAKEN2(hybrid_assemblies.hybrid_contigs)
    
    tb_summary_ch    = channel.of()
    vsnp_summary_ch  = channel.of()
    sistr_summary_ch = channel.of()
    amr_summary_ch   = channel.of()
    mlst_summary_ch  = channel.of()
    vir_summary_ch   = channel.of()
    
    if (params.kraken_db) {
        tb_inputs = hybrid_inputs.join(kraken_hybrid_out.kraken_reports)
                     .filter { id, r1, r2, ont, rep -> decideOrganism(rep.toString()) == 'mbovis' }
                     .map    { id, r1, r2, ont, rep -> tuple(id, r1, r2) }
        TB_PROFILER(tb_inputs)
        tb_mqc = TB_PROFILER.out.tbprofiler_mqc.map { id, f -> f }.collect()
        TB_SUMMARY(tb_mqc)
        tb_summary_ch = TB_SUMMARY.out.tb_summary
        
        if (params.vsnp_ref) {
            VSNP(tb_inputs)
            vsnp_tsv = VSNP.out.vsnp_result.map { id, f -> f }.collect()
            VSNP_SUMMARY(vsnp_tsv)
            vsnp_summary_ch = VSNP_SUMMARY.out.vsnp_summary.collect().flatten()
        }
        
        sal_inputs = hybrid_assemblies.hybrid_contigs.join(kraken_hybrid_out.kraken_reports)
                      .filter { id, c, rep -> decideOrganism(rep.toString()) == 'salmonella' }
                      .map    { id, c, rep -> tuple(id, c) }
        SISTR2(sal_inputs)
        sistr_tsv = SISTR2.out.sistr2_tsv.map { id, f -> f }.collect()
        SISTR_SUMMARY(sistr_tsv)
        sistr_summary_ch = SISTR_SUMMARY.out.sistr_summary
        
        AMRFINDER(sal_inputs)
        amr_reports = AMRFINDER.out.amr_tsv.map { id, f -> f }.collect()
        AMR_SUMMARY(amr_reports)
        amr_summary_ch = AMR_SUMMARY.out.amr_summary
        
        lis_inputs = hybrid_assemblies.hybrid_contigs.join(kraken_hybrid_out.kraken_reports)
                      .filter { id, c, rep -> decideOrganism(rep.toString()) == 'listeria' }
                      .map    { id, c, rep -> tuple(id, c) }
        MLST(lis_inputs)
        mlst_tsv = MLST.out.mlst_tsv.map { id, f -> f }.collect()
        MLST_SUMMARY(mlst_tsv)
        mlst_summary_ch = MLST_SUMMARY.out.mlst_summary
        
        VIRULENCEFINDER(lis_inputs)
        vir_jsons = VIRULENCEFINDER.out.vir_json.map { id, f -> f }.collect()
        VIRULENCE_SUMMARY(vir_jsons)
        vir_summary_ch = VIRULENCE_SUMMARY.out.vir_summary
    }
    
    quast_reports = quast_hybrid_out.quast_reports.map { p -> p }
    
    all_reports = quast_reports
        .mix(tb_summary_ch, vsnp_summary_ch, sistr_summary_ch, amr_summary_ch, mlst_summary_ch, vir_summary_ch)
        .collect()
    
    MULTIQC(all_reports, file("${baseDir}/multiqc_config.yaml"))
    
    log.info "Hybrid assembly pipeline complete. Check results/hybrid_assembly/ for contigs."
}

// ===================== TEST WORKFLOWS =====================
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

workflow test_VIRULENCEFINDER {
    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    trimmed = FASTP(raw_pairs)
    assemblies = SPADES(trimmed.trimmed_reads)
    
    assemblies.contigs
        .map { it.toString() }
        .view()
    
    VIRULENCEFINDER(assemblies.contigs)
}

workflow test_VSNP {
    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    trimmed = FASTP(raw_pairs)
    
    trimmed.trimmed_reads
        .map { it.toString() }
        .view()
    
    VSNP(trimmed.trimmed_reads)
}
