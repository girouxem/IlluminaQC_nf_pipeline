#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ===================== INCLUDE MODULES =====================
include { SANITY_CHECK; FASTQC_RAW; FASTP; FASTQC_TRIMMED } from './modules/qc'
include { RETENTION_TRACK; RETENTION_VIZ } from './modules/retention'
include { SPADES; FILTER_CONTIGS; BUSCO; BUSCO_SUMMARY } from './modules/assembly'

// Import QUAST with original name and aliases for metadata_driven
include { QUAST; QUAST as QUAST_ILL; QUAST as QUAST_ONT; QUAST as QUAST_HYB } from './modules/assembly'

// Import KRAKEN2 with original name and aliases
include { KRAKEN2; KRAKEN2 as KRAKEN2_ILL; KRAKEN2 as KRAKEN2_ONT; KRAKEN2 as KRAKEN2_HYB } from './modules/typing'

// Import TB-Profiler and VSNP with original names and aliases
include { TB_PROFILER; TB_SUMMARY; VSNP; VSNP_SUMMARY } from './modules/typing'
include { TB_PROFILER as TB_PROFILER_ILL; TB_SUMMARY as TB_SUMMARY_ILL; VSNP as VSNP_ILL; VSNP_SUMMARY as VSNP_SUMMARY_ILL } from './modules/typing'
include { TB_PROFILER as TB_PROFILER_ONT; TB_SUMMARY as TB_SUMMARY_ONT; VSNP as VSNP_ONT; VSNP_SUMMARY as VSNP_SUMMARY_ONT } from './modules/typing'
include { TB_PROFILER as TB_PROFILER_HYB; TB_SUMMARY as TB_SUMMARY_HYB; VSNP as VSNP_HYB; VSNP_SUMMARY as VSNP_SUMMARY_HYB } from './modules/typing'

// Import SISTR2, AMRFinder, MLST, VirulenceFinder with original names and aliases
include { SISTR2; SISTR_SUMMARY; AMRFINDER; AMR_SUMMARY; MLST; MLST_SUMMARY; VIRULENCEFINDER; VIRULENCE_SUMMARY } from './modules/typing'
include { SISTR2 as SISTR2_ILL; SISTR_SUMMARY as SISTR_SUMMARY_ILL; AMRFINDER as AMRFINDER_ILL; AMR_SUMMARY as AMR_SUMMARY_ILL; MLST as MLST_ILL; MLST_SUMMARY as MLST_SUMMARY_ILL; VIRULENCEFINDER as VIRULENCEFINDER_ILL; VIRULENCE_SUMMARY as VIRULENCE_SUMMARY_ILL } from './modules/typing'
include { SISTR2 as SISTR2_ONT; SISTR_SUMMARY as SISTR_SUMMARY_ONT; AMRFINDER as AMRFINDER_ONT; AMR_SUMMARY as AMR_SUMMARY_ONT; MLST as MLST_ONT; MLST_SUMMARY as MLST_SUMMARY_ONT; VIRULENCEFINDER as VIRULENCEFINDER_ONT; VIRULENCE_SUMMARY as VIRULENCE_SUMMARY_ONT } from './modules/typing'
include { SISTR2 as SISTR2_HYB; SISTR_SUMMARY as SISTR_SUMMARY_HYB; AMRFINDER as AMRFINDER_HYB; AMR_SUMMARY as AMR_SUMMARY_HYB; MLST as MLST_HYB; MLST_SUMMARY as MLST_SUMMARY_HYB; VIRULENCEFINDER as VIRULENCEFINDER_HYB; VIRULENCE_SUMMARY as VIRULENCE_SUMMARY_HYB } from './modules/typing'

// Import Nanopore processes with original names and aliases for CHECKM
include { DORADO; PORECHOP; FLYE; MEDAKA; CHECKM; CHECKM as CHECKM_ONT; CHECKM as CHECKM_HYB } from './modules/nanopore'
include { UNICYCLER } from './modules/hybrid'
include { CHECKPOINT_REPORT; ASSEMBLY_VIZ; READ_LEN_HIST } from './modules/viz'
include { MULTIQC } from './modules/multiqc'
include { LOAD_METADATA } from './modules/metadata'

// ===================== PARAMETERS =====================
params.test = params.test ?: null
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
params.samples_csv    = params.samples_csv    ?: null

// ===================== HELPER FUNCTIONS =====================
def decideOrganism(String kraken_report_path) {
    def lines = new File(kraken_report_path).readLines()
    def species = lines.find { line -> line.contains("\tS\t") } ?: ''
    if (species.contains("Mycobacterium")) return "mbovis"
    if (species.contains("Salmonella"))    return "salmonella"
    if (species.contains("Listeria"))      return "listeria"
    if (species.toLowerCase().contains("virus")) return "virus"
    return "bacteria"
}

// ===================== MAIN WORKFLOW ROUTER =====================
workflow {
    if (params.test == "metadata") {
        test_LOAD_METADATA()
    } else if (params.test == "sanity_check") {
        test_SANITY_CHECK()
    } else if (params.test == "fastqc_raw") {
        test_FASTQC_RAW()
    } else if (params.test == "fastp") {
        test_FASTP()
    } else if (params.test == "fastqc_trimmed") {
        test_FASTQC_TRIMMED()
    } else if (params.test == "retention") {
        test_RETENTION_TRACK()
    } else if (params.test == "spades") {
        test_SPADES()
    } else if (params.test == "quast") {
        test_QUAST()
    } else if (params.test == "busco") {
        test_BUSCO()
    } else if (params.test == "kraken2") {
        test_KRAKEN2()
    } else if (params.test == "tbprofiler") {
        test_TB_PROFILER()
    } else if (params.test == "kraken2_routing") {
        test_KRAKEN2_ROUTING()
    } else if (params.test == "virulencefinder") {
        test_VIRULENCEFINDER()
    } else if (params.test == "vsnp") {
        test_VSNP()
    } else if (params.run_nanopore) {
        nanopore_pipeline()
    } else if (params.run_hybrid) {
        hybrid_assembly()
    } else if (params.samples_csv) {
        metadata_driven()
    } else {
        // Standard Illumina Pipeline
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
                         .filter { id, _r1, _r2, rep -> decideOrganism(rep.toString()) == 'mbovis' }
                         .map    { id, r1, r2, _rep -> tuple(id, r1, r2) }
            TB_PROFILER(tb_inputs)
            tb_mqc_reports = TB_PROFILER.out.tbprofiler_mqc.map { _id, f -> f }.collect().unique { file -> file.name }
            TB_SUMMARY(tb_mqc_reports)

            if (params.vsnp_ref) {
                VSNP(tb_inputs)
                vsnp_tsv = VSNP.out.vsnp_result.map { _id, f -> f }.collect().unique { file -> file.name }
                VSNP_SUMMARY(vsnp_tsv)
                vsnp_summary_file = VSNP_SUMMARY.out.vsnp_summary.collect().flatten()
            }

            sal_inputs = assemblies.contigs.join(kraken_out.kraken_reports)
                          .filter { _id, _c, rep -> decideOrganism(rep.toString()) == 'salmonella' }
                          .map    { id, c, _rep -> tuple(id, c) }
            SISTR2(sal_inputs)
            sistr_tsv = SISTR2.out.sistr2_tsv.map { _id, f -> f }.collect().unique { file -> file.name }
            SISTR_SUMMARY(sistr_tsv)

            AMRFINDER(sal_inputs)
            amr_reports = AMRFINDER.out.amr_tsv.map { _id, f -> f }.collect().unique { file -> file.name }
            AMR_SUMMARY(amr_reports)

            lis_inputs = assemblies.contigs.join(kraken_out.kraken_reports)
                          .filter { _id, _c, rep -> decideOrganism(rep.toString()) == 'listeria' }
                          .map    { id, c, _rep -> tuple(id, c) }
            MLST(lis_inputs)
            mlst_tsv = MLST.out.mlst_tsv.map { _id, f -> f }.collect().unique { file -> file.name }
            MLST_SUMMARY(mlst_tsv)

            VIRULENCEFINDER(lis_inputs)
            vir_jsons = VIRULENCEFINDER.out.vir_json.map { _id, f -> f }.collect().unique { file -> file.name }
            VIRULENCE_SUMMARY(vir_jsons)
        }

        raw_reports     = raw_fastqc.fastqc_zip.map { _id, f -> f }
        trimmed_reports = trimmed.fastp_json.map { _id, f -> f }
        post_reports    = trimmed_fastqc.fastqc_zip_trim.map { _id, f -> f }
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

    needs_basecalling = ont_inputs.filter { _id, f -> f.name.endsWith('.pod5') || f.name.endsWith('.fast5') }
    already_basecalled = ont_inputs.filter { _id, f -> f.name.endsWith('.fastq') || f.name.endsWith('.fastq.gz') }
    
    basecalled = DORADO(needs_basecalling)
    all_ont_reads = basecalled.basecalled_fastq.mix(already_basecalled)
    
    trimmed_ont = PORECHOP(all_ont_reads)
    ont_assemblies = FLYE(trimmed_ont.trimmed_fastq)
    
    polished_inputs = ont_assemblies.contigs.join(trimmed_ont.trimmed_fastq)
    polished_ont = MEDAKA(polished_inputs)
    
    quast_ont_out = QUAST(polished_ont.polished_contigs)
    _checkm_ont_out = CHECKM(polished_ont.polished_contigs)
    
    kraken_ont_out = KRAKEN2(polished_ont.polished_contigs)
    
    tb_summary_ch    = channel.of()
    vsnp_summary_ch  = channel.of()
    sistr_summary_ch = channel.of()
    amr_summary_ch   = channel.of()
    mlst_summary_ch  = channel.of()
    vir_summary_ch   = channel.of()
    
    if (params.kraken_db) {
        tb_ont_inputs = trimmed_ont.trimmed_fastq.join(kraken_ont_out.kraken_reports)
                         .filter { _id, _fq, rep -> decideOrganism(rep.toString()) == 'mbovis' }
                         .map    { id, fq, _rep -> tuple(id, fq, fq) }
        TB_PROFILER(tb_ont_inputs)
        tb_mqc_reports = TB_PROFILER.out.tbprofiler_mqc.map { _id, f -> f }.collect().unique { file -> file.name }
        TB_SUMMARY(tb_mqc_reports)
        tb_summary_ch = TB_SUMMARY.out.tb_summary
        
        if (params.vsnp_ref) {
            VSNP(tb_ont_inputs)
            vsnp_tsv = VSNP.out.vsnp_result.map { _id, f -> f }.collect().unique { file -> file.name }
            VSNP_SUMMARY(vsnp_tsv)
            vsnp_summary_ch = VSNP_SUMMARY.out.vsnp_summary.collect().flatten()
        }
        
        sal_ont_inputs = polished_ont.polished_contigs.join(kraken_ont_out.kraken_reports)
                          .filter { _id, _c, rep -> decideOrganism(rep.toString()) == 'salmonella' }
                          .map    { id, c, _rep -> tuple(id, c) }
        SISTR2(sal_ont_inputs)
        sistr_tsv = SISTR2.out.sistr2_tsv.map { _id, f -> f }.collect().unique { file -> file.name }
        SISTR_SUMMARY(sistr_tsv)
        sistr_summary_ch = SISTR_SUMMARY.out.sistr_summary
        
        AMRFINDER(sal_ont_inputs)
        amr_reports = AMRFINDER.out.amr_tsv.map { _id, f -> f }.collect().unique { file -> file.name }
        AMR_SUMMARY(amr_reports)
        amr_summary_ch = AMR_SUMMARY.out.amr_summary
        
        lis_ont_inputs = polished_ont.polished_contigs.join(kraken_ont_out.kraken_reports)
                          .filter { _id, _c, rep -> decideOrganism(rep.toString()) == 'listeria' }
                          .map    { id, c, _rep -> tuple(id, c) }
        MLST(lis_ont_inputs)
        mlst_tsv = MLST.out.mlst_tsv.map { _id, f -> f }.collect().unique { file -> file.name }
        MLST_SUMMARY(mlst_tsv)
        mlst_summary_ch = MLST_SUMMARY.out.mlst_summary
        
        VIRULENCEFINDER(lis_ont_inputs)
        vir_jsons = VIRULENCEFINDER.out.vir_json.map { _id, f -> f }.collect().unique { file -> file.name }
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

    needs_basecalling = ont_inputs.filter { _id, f -> f.name.endsWith('.pod5') || f.name.endsWith('.fast5') }
    already_basecalled = ont_inputs.filter { _id, f -> f.name.endsWith('.fastq') || f.name.endsWith('.fastq.gz') }
    basecalled = DORADO(needs_basecalling)
    all_ont_reads = basecalled.basecalled_fastq.mix(already_basecalled)
    ont_trimmed = PORECHOP(all_ont_reads)
    
    hybrid_inputs = illumina_trimmed.trimmed_reads.join(ont_trimmed.trimmed_fastq)
    
    hybrid_inputs.map { id, _r1, _r2, _ont -> "HYBRID: ${id} has Illumina + Nanopore" }.view()
    
    hybrid_assemblies = UNICYCLER(hybrid_inputs)
    
    quast_hybrid_out = QUAST(hybrid_assemblies.hybrid_contigs)
    _checkm_hybrid_out = CHECKM(hybrid_assemblies.hybrid_contigs)
    
    kraken_hybrid_out = KRAKEN2(hybrid_assemblies.hybrid_contigs)
    
    tb_summary_ch    = channel.of()
    vsnp_summary_ch  = channel.of()
    sistr_summary_ch = channel.of()
    amr_summary_ch   = channel.of()
    mlst_summary_ch  = channel.of()
    vir_summary_ch   = channel.of()
    
    if (params.kraken_db) {
        tb_inputs = hybrid_inputs.join(kraken_hybrid_out.kraken_reports)
                     .filter { _id, _r1, _r2, _ont, rep -> decideOrganism(rep.toString()) == 'mbovis' }
                     .map    { id, r1, r2, _ont, _rep -> tuple(id, r1, r2) }
        TB_PROFILER(tb_inputs)
        tb_mqc = TB_PROFILER.out.tbprofiler_mqc.map { _id, f -> f }.collect().unique { file -> file.name }
        TB_SUMMARY(tb_mqc)
        tb_summary_ch = TB_SUMMARY.out.tb_summary
        
        if (params.vsnp_ref) {
            VSNP(tb_inputs)
            vsnp_tsv = VSNP.out.vsnp_result.map { _id, f -> f }.collect().unique { file -> file.name }
            VSNP_SUMMARY(vsnp_tsv)
            vsnp_summary_ch = VSNP_SUMMARY.out.vsnp_summary.collect().flatten()
        }
        
        sal_inputs = hybrid_assemblies.hybrid_contigs.join(kraken_hybrid_out.kraken_reports)
                      .filter { _id, _c, rep -> decideOrganism(rep.toString()) == 'salmonella' }
                      .map    { id, c, _rep -> tuple(id, c) }
        SISTR2(sal_inputs)
        sistr_tsv = SISTR2.out.sistr2_tsv.map { _id, f -> f }.collect().unique { file -> file.name }
        SISTR_SUMMARY(sistr_tsv)
        sistr_summary_ch = SISTR_SUMMARY.out.sistr_summary
        
        AMRFINDER(sal_inputs)
        amr_reports = AMRFINDER.out.amr_tsv.map { _id, f -> f }.collect().unique { file -> file.name }
        AMR_SUMMARY(amr_reports)
        amr_summary_ch = AMR_SUMMARY.out.amr_summary
        
        lis_inputs = hybrid_assemblies.hybrid_contigs.join(kraken_hybrid_out.kraken_reports)
                      .filter { _id, _c, rep -> decideOrganism(rep.toString()) == 'listeria' }
                      .map    { id, c, _rep -> tuple(id, c) }
        MLST(lis_inputs)
        mlst_tsv = MLST.out.mlst_tsv.map { _id, f -> f }.collect().unique { file -> file.name }
        MLST_SUMMARY(mlst_tsv)
        mlst_summary_ch = MLST_SUMMARY.out.mlst_summary
        
        VIRULENCEFINDER(lis_inputs)
        vir_jsons = VIRULENCEFINDER.out.vir_json.map { _id, f -> f }.collect().unique { file -> file.name }
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

// ===================== METADATA DRIVEN WORKFLOW =====================
workflow metadata_driven {
    log.info "=== Running Metadata-Driven Pipeline ==="
    log.info "Samples CSV: ${params.samples_csv}"
    
    // 1. Load metadata and capture the output channels
    metadata_out = LOAD_METADATA(file(params.samples_csv))
    
    // 2. Parse the channels TSV directly from the process output channel
    sample_channels = metadata_out.channels_tsv
        .splitCsv(sep: '\t', header: true)
        .map { row -> tuple(row.sample_id, row.data_type, row.organism, 
                           row.illumina_r1, row.illumina_r2, row.nanopore) }
    
    // 3. Route based on data_type
    illumina_samples = sample_channels.filter { _id, type, _org, _r1, _r2, _ont -> 
        type == "illumina" && _r1 && _r2 
    }
    nanopore_samples = sample_channels.filter { _id, type, _org, _r1, _r2, _ont -> 
        type == "nanopore" && _ont 
    }
    hybrid_samples = sample_channels.filter { _id, type, _org, _r1, _r2, _ont -> 
        type == "hybrid" 
    }

    // ---- ILLUMINA BRANCH ----
    illumina_input = illumina_samples.map { id, _type, _org, r1, r2, _ont -> 
        tuple(id, [file(r1), file(r2)]) 
    }
    illumina_trimmed = FASTP(illumina_input)

    illumina_assemblies = SPADES(illumina_trimmed.trimmed_reads)
    illumina_filtered   = FILTER_CONTIGS(illumina_assemblies.contigs)
    illumina_quast_out  = QUAST_ILL(illumina_assemblies.contigs)
    illumina_busco_out  = BUSCO(illumina_filtered.filtered_contigs)
    illumina_busco_collected = illumina_busco_out.busco_reports.collect()
    BUSCO_SUMMARY(illumina_busco_collected)

    // ---- NANOPORE BRANCH ----
    ont_input = nanopore_samples.map { id, _type, _org, _r1, _r2, ont -> 
        tuple(id, file(ont)) 
    }
    ont_trimmed = PORECHOP(ont_input)
    ont_assemblies = FLYE(ont_trimmed.trimmed_fastq)
    
    ont_polished_inputs = ont_assemblies.contigs.join(ont_trimmed.trimmed_fastq)
    ont_polished = MEDAKA(ont_polished_inputs)
    
    ont_quast_out = QUAST_ONT(ont_polished.polished_contigs)
    ont_checkm_out = CHECKM_ONT(ont_polished.polished_contigs)

    // ---- HYBRID BRANCH ----
    hybrid_input = hybrid_samples.map { id, _type, _org, r1, r2, ont -> 
        tuple(id, file(r1), file(r2), file(ont)) 
    }
    hybrid_assemblies = UNICYCLER(hybrid_input)
    
    hybrid_quast_out = QUAST_HYB(hybrid_assemblies.hybrid_contigs)
    hybrid_checkm_out = CHECKM_HYB(hybrid_assemblies.hybrid_contigs)

    // ---- Kraken2 classification for all assemblies ----
    kraken_illumina = KRAKEN2_ILL(illumina_assemblies.contigs)
    kraken_ont      = KRAKEN2_ONT(ont_polished.polished_contigs)
    kraken_hybrid   = KRAKEN2_HYB(hybrid_assemblies.hybrid_contigs)

    // ---- Initialize empty channels for per-sample typing files ----
    tb_ill_mqc       = channel.of()
    vsnp_ill_tsv     = channel.of()
    sistr_ill_tsv    = channel.of()
    amr_ill_reports  = channel.of()
    mlst_ill_tsv     = channel.of()
    vir_ill_jsons    = channel.of()

    tb_ont_mqc       = channel.of()
    vsnp_ont_tsv     = channel.of()
    sistr_ont_tsv    = channel.of()
    amr_ont_reports  = channel.of()
    mlst_ont_tsv     = channel.of()
    vir_ont_jsons    = channel.of()

    tb_hyb_mqc       = channel.of()
    vsnp_hyb_tsv     = channel.of()
    sistr_hyb_tsv    = channel.of()
    amr_hyb_reports  = channel.of()
    mlst_hyb_tsv     = channel.of()
    vir_hyb_jsons    = channel.of()

    // ---- Route Illumina assemblies to typing tools ----
    if (params.kraken_db) {
        tb_ill_inputs = illumina_trimmed.trimmed_reads.join(kraken_illumina.kraken_reports)
                         .filter { _id, _r1, _r2, rep -> decideOrganism(rep.toString()) == 'mbovis' }
                         .map    { id, r1, r2, _rep -> tuple(id, r1, r2) }
        TB_PROFILER_ILL(tb_ill_inputs)
        tb_ill_mqc = TB_PROFILER_ILL.out.tbprofiler_mqc.map { _id, f -> f }.collect()
        
        if (params.vsnp_ref) {
            VSNP_ILL(tb_ill_inputs)
            vsnp_ill_tsv = VSNP_ILL.out.vsnp_result.map { _id, f -> f }.collect()
        }

        sal_ill_inputs = illumina_assemblies.contigs.join(kraken_illumina.kraken_reports)
                          .filter { _id, _c, rep -> decideOrganism(rep.toString()) == 'salmonella' }
                          .map    { id, c, _rep -> tuple(id, c) }
        SISTR2_ILL(sal_ill_inputs)
        sistr_ill_tsv = SISTR2_ILL.out.sistr2_tsv.map { _id, f -> f }.collect()
        
        AMRFINDER_ILL(sal_ill_inputs)
        amr_ill_reports = AMRFINDER_ILL.out.amr_tsv.map { _id, f -> f }.collect()

        lis_ill_inputs = illumina_assemblies.contigs.join(kraken_illumina.kraken_reports)
                          .filter { _id, _c, rep -> decideOrganism(rep.toString()) == 'listeria' }
                          .map    { id, c, _rep -> tuple(id, c) }
        MLST_ILL(lis_ill_inputs)
        mlst_ill_tsv = MLST_ILL.out.mlst_tsv.map { _id, f -> f }.collect()
        
        VIRULENCEFINDER_ILL(lis_ill_inputs)
        vir_ill_jsons = VIRULENCEFINDER_ILL.out.vir_json.map { _id, f -> f }.collect()
    }

    // ---- Route ONT assemblies to typing tools ----
    if (params.kraken_db) {
        tb_ont_inputs = ont_trimmed.trimmed_fastq.join(kraken_ont.kraken_reports)
                         .filter { _id, _fq, rep -> decideOrganism(rep.toString()) == 'mbovis' }
                         .map    { id, fq, _rep -> tuple(id, fq, fq) }
        TB_PROFILER_ONT(tb_ont_inputs)
        tb_ont_mqc = TB_PROFILER_ONT.out.tbprofiler_mqc.map { _id, f -> f }.collect()
        
        if (params.vsnp_ref) {
            VSNP_ONT(tb_ont_inputs)
            vsnp_ont_tsv = VSNP_ONT.out.vsnp_result.map { _id, f -> f }.collect()
        }

        sal_ont_inputs = ont_polished.polished_contigs.join(kraken_ont.kraken_reports)
                          .filter { _id, _c, rep -> decideOrganism(rep.toString()) == 'salmonella' }
                          .map    { id, c, _rep -> tuple(id, c) }
        SISTR2_ONT(sal_ont_inputs)
        sistr_ont_tsv = SISTR2_ONT.out.sistr2_tsv.map { _id, f -> f }.collect()
        
        AMRFINDER_ONT(sal_ont_inputs)
        amr_ont_reports = AMRFINDER_ONT.out.amr_tsv.map { _id, f -> f }.collect()

        lis_ont_inputs = ont_polished.polished_contigs.join(kraken_ont.kraken_reports)
                          .filter { _id, _c, rep -> decideOrganism(rep.toString()) == 'listeria' }
                          .map    { id, c, _rep -> tuple(id, c) }
        MLST_ONT(lis_ont_inputs)
        mlst_ont_tsv = MLST_ONT.out.mlst_tsv.map { _id, f -> f }.collect()
        
        VIRULENCEFINDER_ONT(lis_ont_inputs)
        vir_ont_jsons = VIRULENCEFINDER_ONT.out.vir_json.map { _id, f -> f }.collect()
    }

    // ---- Route Hybrid assemblies to typing tools ----
    if (params.kraken_db) {
        tb_hyb_inputs = hybrid_input.join(kraken_hybrid.kraken_reports)
                         .filter { _id, _r1, _r2, _ont, rep -> decideOrganism(rep.toString()) == 'mbovis' }
                         .map    { id, r1, r2, _ont, _rep -> tuple(id, r1, r2) }
        TB_PROFILER_HYB(tb_hyb_inputs)
        tb_hyb_mqc = TB_PROFILER_HYB.out.tbprofiler_mqc.map { _id, f -> f }.collect()
        
        if (params.vsnp_ref) {
            VSNP_HYB(tb_hyb_inputs)
            vsnp_hyb_tsv = VSNP_HYB.out.vsnp_result.map { _id, f -> f }.collect()
        }

        sal_hyb_inputs = hybrid_assemblies.hybrid_contigs.join(kraken_hybrid.kraken_reports)
                          .filter { _id, _c, rep -> decideOrganism(rep.toString()) == 'salmonella' }
                          .map    { id, c, _rep -> tuple(id, c) }
        SISTR2_HYB(sal_hyb_inputs)
        sistr_hyb_tsv = SISTR2_HYB.out.sistr2_tsv.map { _id, f -> f }.collect()
        
        AMRFINDER_HYB(sal_hyb_inputs)
        amr_hyb_reports = AMRFINDER_HYB.out.amr_tsv.map { _id, f -> f }.collect()

        lis_hyb_inputs = hybrid_assemblies.hybrid_contigs.join(kraken_hybrid.kraken_reports)
                          .filter { _id, _c, rep -> decideOrganism(rep.toString()) == 'listeria' }
                          .map    { id, c, _rep -> tuple(id, c) }
        MLST_HYB(lis_hyb_inputs)
        mlst_hyb_tsv = MLST_HYB.out.mlst_tsv.map { _id, f -> f }.collect()
        
        VIRULENCEFINDER_HYB(lis_hyb_inputs)
        vir_hyb_jsons = VIRULENCEFINDER_HYB.out.vir_json.map { _id, f -> f }.collect()
    }

    // ---- Combine all branches and run summaries ONCE ----
    // Initialize empty channels for the final summary files
    tb_summary_ch    = channel.of()
    vsnp_summary_ch  = channel.of()
    sistr_summary_ch = channel.of()
    amr_summary_ch   = channel.of()
    mlst_summary_ch  = channel.of()
    vir_summary_ch   = channel.of()

    if (params.kraken_db) {
        // TB
        TB_SUMMARY(tb_ill_mqc.mix(tb_ont_mqc, tb_hyb_mqc).collect())
        tb_summary_ch = TB_SUMMARY.out.tb_summary.collect().flatten()
        
        // VSNP
        if (params.vsnp_ref) {
            VSNP_SUMMARY(vsnp_ill_tsv.mix(vsnp_ont_tsv, vsnp_hyb_tsv).collect())
            vsnp_summary_ch = VSNP_SUMMARY.out.vsnp_summary.collect().flatten()
        }
        
        // SISTR
        SISTR_SUMMARY(sistr_ill_tsv.mix(sistr_ont_tsv, sistr_hyb_tsv).collect())
        sistr_summary_ch = SISTR_SUMMARY.out.sistr_summary.collect().flatten()
        
        // AMR
        AMR_SUMMARY(amr_ill_reports.mix(amr_ont_reports, amr_hyb_reports).collect())
        amr_summary_ch = AMR_SUMMARY.out.amr_summary.collect().flatten()
        
        // MLST
        MLST_SUMMARY(mlst_ill_tsv.mix(mlst_ont_tsv, mlst_hyb_tsv).collect())
        mlst_summary_ch = MLST_SUMMARY.out.mlst_summary.collect().flatten()
        
        // Virulence
        VIRULENCE_SUMMARY(vir_ill_jsons.mix(vir_ont_jsons, vir_hyb_jsons).collect())
        vir_summary_ch = VIRULENCE_SUMMARY.out.vir_summary.collect().flatten()
    }

    // ---- Collect QC reports for MultiQC ----
    fastp_reports = illumina_trimmed.fastp_json.map { _id, f -> f }.collect().ifEmpty([])
    
    illumina_quast_reports = illumina_quast_out.quast_reports.map { p -> p }.collect().ifEmpty([])
    ont_quast_reports      = ont_quast_out.quast_reports.map { p -> p }.collect().ifEmpty([])
    hybrid_quast_reports   = hybrid_quast_out.quast_reports.map { p -> p }.collect().ifEmpty([])
    
    busco_summary_file = BUSCO_SUMMARY.out.busco_summary.collect().flatten().ifEmpty([])

    ont_checkm_reports    = ont_checkm_out.checkm_result.map { _id, f -> f }.collect().ifEmpty([])
    hybrid_checkm_reports = hybrid_checkm_out.checkm_result.map { _id, f -> f }.collect().ifEmpty([])

    // ---- Feed everything to MultiQC ----
    all_reports = metadata_out.metadata_tsv
        .mix(fastp_reports)
        .mix(illumina_quast_reports, ont_quast_reports, hybrid_quast_reports)
        .mix(busco_summary_file)
        .mix(ont_checkm_reports, hybrid_checkm_reports)
        .mix(tb_summary_ch, vsnp_summary_ch, sistr_summary_ch, amr_summary_ch, mlst_summary_ch, vir_summary_ch)
        .collect()
    
    MULTIQC(all_reports, file("${baseDir}/multiqc_config.yaml"))
    
    log.info "Metadata-driven pipeline complete."
}


// ===================== TEST WORKFLOWS =====================
workflow test_LOAD_METADATA {
    log.info "=== Testing Metadata CSV Parsing ==="
    log.info "CSV file: ${params.samples_csv}"
    
    if (!file(params.samples_csv).exists()) {
        log.error "ERROR: CSV file not found: ${params.samples_csv}"
        System.exit(1)
    }
    
    metadata_out = LOAD_METADATA(file(params.samples_csv))
    
    metadata_out.metadata_tsv.view { f ->
        "=== METADATA TSV ===\n" + f.text
    }
    
    metadata_out.channels_tsv.view { f ->
        "=== CHANNELS TSV ===\n" + f.text
    }
    
    log.info "=== Metadata test complete ==="
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
    routed.map { id, _contigs, rep ->
        def org = decideOrganism(rep.toString())
        "ROUTING: ${id} -> ${org}"
    }.view()

    tb_inputs = trimmed.trimmed_reads.join(kraken_out.kraken_reports)
                 .filter { _id, _r1, _r2, rep -> decideOrganism(rep.toString()) == 'mbovis' }
                 .map    { id, r1, r2, _rep -> tuple(id, r1, r2) }
    tb_inputs.map { id, _r1, _r2 -> "TB_PROFILER will run on: ${id}" }.view()
    TB_PROFILER(tb_inputs)

    sal_inputs = assemblies.contigs.join(kraken_out.kraken_reports)
                  .filter { _id, _c, rep -> decideOrganism(rep.toString()) == 'salmonella' }
                  .map    { id, c, _rep -> tuple(id, c) }
    sal_inputs.map { id, _c -> "SISTR2 + AMRFINDER will run on: ${id}" }.view()
    SISTR2(sal_inputs)
    AMRFINDER(sal_inputs)

    lis_inputs = assemblies.contigs.join(kraken_out.kraken_reports)
                  .filter { _id, _c, rep -> decideOrganism(rep.toString()) == 'listeria' }
                  .map    { id, c, _rep -> tuple(id, c) }
    lis_inputs.map { id, _c -> "MLST + VIRULENCEFINDER will run on: ${id}" }.view()
    MLST(lis_inputs)
    VIRULENCEFINDER(lis_inputs)
}

workflow test_VIRULENCEFINDER {
    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    trimmed = FASTP(raw_pairs)
    assemblies = SPADES(trimmed.trimmed_reads)
    
    assemblies.contigs
        .map { val -> val.toString() }
        .view()
    
    VIRULENCEFINDER(assemblies.contigs)
}

workflow test_VSNP {
    raw_pairs = channel.fromFilePairs(params.reads, flat: false)
    trimmed = FASTP(raw_pairs)
    
    trimmed.trimmed_reads
        .map { val -> val.toString() }
        .view()
    
    VSNP(trimmed.trimmed_reads)
}
