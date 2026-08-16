// modules/retention.nf
// Contains: RETENTION_TRACK, RETENTION_VIZ

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
