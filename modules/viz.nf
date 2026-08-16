// modules/viz.nf
// Contains: CHECKPOINT_REPORT, ASSEMBLY_VIZ, READ_LEN_HIST

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
